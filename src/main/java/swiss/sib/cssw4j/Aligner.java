package swiss.sib.cssw4j;

import static jdk.incubator.vector.VectorOperators.GT;
import static jdk.incubator.vector.VectorOperators.MAX;
import static jdk.incubator.vector.VectorOperators.SADD;
import static jdk.incubator.vector.VectorOperators.SSUB;
import static jdk.incubator.vector.VectorOperators.SUSUB;

import java.util.Arrays;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * Striped Smith-Waterman aligner, a port of {@code ssw_align} from the Complete Striped Smith-Waterman library
 * (Zhao et al. 2013, Farrar 2007) to the Java Vector API.
 * <p>
 * An aligner owns all the scratch memory the kernels need and reuses it between calls, so steady state alignment does
 * not allocate beyond the result. It is therefore <b>not thread safe</b>: use one instance per thread. {@link Profile}
 * and {@link ScoringMatrix} can be shared.
 * <p>
 * As in the C library, scores are exact Smith-Waterman scores only for {@code gapOpen > gapExtend}: with
 * {@code gapOpen <= gapExtend} the lazy-F loop may stop propagating gaps along the query too early.
 */
public final class Aligner {
	/** Always compute begin positions and the cigar (bit 8 of the C flag). */
	public static final int FLAG_ALWAYS = 0x01;
	/** Compute the cigar when {@code score1 >= filters} (bit 7). */
	public static final int FLAG_SCORE_FILTER = 0x02;
	/** Compute the cigar when both aligned spans are {@code <= filterd} (bit 6). */
	public static final int FLAG_DISTANCE_FILTER = 0x04;
	/** Compute the begin positions (bit 5). */
	public static final int FLAG_BEGIN = 0x08;

	private static final VectorSpecies<Byte> BS = Lanes.BYTE;
	private static final VectorSpecies<Short> WS = Lanes.SHORT;
	private static final int NEG_INF = Integer.MIN_VALUE / 2;
	/*
	 * Banded traceback directions, packed in one byte per cell. The low 3 bits hold the H direction with the C codes
	 * (1: diagonal, 2/3: from E extended/opened, 4/5: from F extended/opened); the two flags record whether E and F
	 * were opened (C codes 3 and 5) or extended (2 and 4).
	 */
	private static final int DIR_E_OPEN = 8, DIR_F_OPEN = 16;
	/** Score 0 of the 8 bit kernel in excess-128 form. */
	private static final byte ZERO = Byte.MIN_VALUE;

	// Striped kernel scratch.
	private byte[] bH1 = new byte[0], bH2 = new byte[0], bE = new byte[0], bHmax = new byte[0];
	private short[] wH1 = new short[0], wH2 = new short[0], wE = new short[0], wHmax = new short[0];
	private int[] maxColumn = new int[0];
	private final byte[] bGaps = new byte[2 * BS.length()];
	private final short[] wGaps = new short[2 * WS.length()];
	private byte[] readReverse = new byte[0];
	private byte[] reverseProfileByte;
	private short[] reverseProfileWord;
	// Results of the last kernel run (the alignment_end pair of the C code).
	private int bestScore, bestRef, bestRead, secondScore, secondRef;
	// Banded alignment scratch.
	private int[] hB = new int[8], eB = new int[8], hC = new int[8];
	private byte[] direction = new byte[1024];
	private int[] cigarBuf = new int[16];

	/**
	 * Full alignment with begin positions and cigar, using {@code readLen / 2} as mask length like {@code ssw_test}.
	 *
	 * @param ref encoded reference
	 */
	public Alignment align(Profile profile, byte[] ref, int gapOpen, int gapExtend) {
		return align(profile, ref, ref.length, gapOpen, gapExtend, FLAG_ALWAYS, 0, 0, profile.readLen / 2);
	}

	/**
	 * Port of {@code ssw_align}.
	 *
	 * @param profile   query profile
	 * @param ref       encoded reference
	 * @param refLen    number of reference letters to use
	 * @param gapOpen   absolute gap open penalty (charged for the first letter of a gap)
	 * @param gapExtend absolute gap extension penalty
	 * @param flag      combination of the {@code FLAG_*} constants; 0 only computes scores and end positions
	 * @param filters   score filter used with {@link #FLAG_SCORE_FILTER}
	 * @param filterd   distance filter used with {@link #FLAG_DISTANCE_FILTER}
	 * @param maskLen   distance between the best and sub-optimal end positions; must be &gt;= 15 to get a
	 *                  sub-optimal score
	 * @throws IllegalStateException when the score overflows the 8 bit kernel and the profile has no 16 bit part
	 */
	public Alignment align(Profile profile, byte[] ref, int refLen, int gapOpen, int gapExtend, int flag, int filters,
			int filterd, int maskLen) {
		if (gapOpen < 0 || gapOpen > 255 || gapExtend < 0 || gapExtend > 255) {
			throw new IllegalArgumentException("gap penalties must be in [0, 255]");
		}
		if (refLen < 0 || refLen > ref.length) {
			throw new IllegalArgumentException("invalid reference length " + refLen);
		}
		final int readLen = profile.readLen;
		boolean word = false;
		// The 8 bit kernel works on signed bytes, larger gap penalties go to the 16 bit kernel.
		final boolean byteGaps = gapOpen <= 127 && gapExtend <= 127;
		if (profile.profileByte != null && (byteGaps || profile.profileWord == null)) {
			if (!byteGaps) {
				throw new IllegalArgumentException("gap penalties above 127 need a profile with ScoreSize.WORD or BOTH");
			}
			swByte(ref, false, refLen, readLen, gapOpen, gapExtend, profile.profileByte, 255, profile.bias, maskLen);
			if (bestScore == 255) {
				if (profile.profileWord == null) {
					throw new IllegalStateException("score overflows 8 bits, build the profile with ScoreSize.BOTH");
				}
				swWord(ref, false, refLen, readLen, gapOpen, gapExtend, profile.profileWord, 0xffff, maskLen);
				word = true;
			}
		} else {
			swWord(ref, false, refLen, readLen, gapOpen, gapExtend, profile.profileWord, 0xffff, maskLen);
			word = true;
		}
		if (bestScore <= 0) {
			return new Alignment(0, 0, -1, 0, -1, 0, 0, null, 0);
		}
		final int score1 = bestScore, refEnd1 = bestRef, readEnd1 = bestRead;
		final int score2 = maskLen >= 15 ? secondScore : 0;
		final int refEnd2 = maskLen >= 15 ? secondRef : -1;
		if (flag == 0 || (flag == FLAG_SCORE_FILTER && score1 < filters)) {
			return new Alignment(score1, score2, -1, refEnd1, -1, readEnd1, refEnd2, null, 0);
		}

		// Find the beginning of the best alignment by aligning the reversed query prefix backwards.
		final int revLen = readEnd1 + 1;
		if (readReverse.length < revLen) {
			readReverse = new byte[Math.max(revLen, readReverse.length * 2)];
		}
		for (int i = 0; i < revLen; i++) {
			readReverse[i] = profile.read[readEnd1 - i];
		}
		final ScoringMatrix m = profile.matrix;
		if (!word) {
			reverseProfileByte = Profile.buildByte(readReverse, revLen, m.matrix(), m.size(), reverseProfileByte);
			swByte(ref, true, refEnd1 + 1, revLen, gapOpen, gapExtend, reverseProfileByte, score1 & 0xff, profile.bias,
					maskLen);
		} else {
			reverseProfileWord = Profile.buildWord(readReverse, revLen, m.matrix(), m.size(), reverseProfileWord);
			swWord(ref, true, refEnd1 + 1, revLen, gapOpen, gapExtend, reverseProfileWord, score1, maskLen);
		}
		final int refBegin1 = bestRef;
		final int readBegin1 = readEnd1 - bestRead;
		int resultFlag = score1 > bestScore ? 2 : 0;

		if ((7 & flag) == 0 || ((2 & flag) != 0 && score1 < filters)
				|| ((4 & flag) != 0 && (refEnd1 - refBegin1 > filterd || readEnd1 - readBegin1 > filterd))) {
			return new Alignment(score1, score2, refBegin1, refEnd1, readBegin1, readEnd1, refEnd2, null, resultFlag);
		}

		// Generate the cigar.
		final int alnRefLen = refEnd1 - refBegin1 + 1;
		final int alnReadLen = readEnd1 - readBegin1 + 1;
		int bandWidth = Math.abs(alnRefLen - alnReadLen) + 1;
		final int fullBand = Math.max(alnRefLen, alnReadLen);
		int[] path;
		while (true) {
			path = bandedSw(ref, refBegin1, profile.read, readBegin1, alnRefLen, alnReadLen, score1, gapOpen,
					gapExtend, bandWidth, m.transposed(), m.size());
			if (path == null) {
				break;
			}
			int cigarScore = cigarScore(path, ref, refBegin1, alnRefLen, profile.read, readBegin1, alnReadLen,
					m.matrix(), m.size(), gapOpen, gapExtend);
			if (cigarScore == score1) {
				break;
			}
			if (bandWidth >= fullBand) {
				path = null;
				break;
			}
			bandWidth = fullBand;
		}
		if (path == null) {
			resultFlag = 1;
		}
		return new Alignment(score1, score2, refBegin1, refEnd1, readBegin1, readEnd1, refEnd2, path, resultFlag);
	}

	/**
	 * 8 bit striped Smith-Waterman (Farrar) with lazy-F loop. Records the best score with its end positions and the
	 * sub-optimal score from the per column maxima.
	 * <p>
	 * The C code keeps scores as unsigned bytes and adds a biased profile. Here scores {@code u} in {@code [0, 255]} are
	 * stored in excess-128 form ({@code u - 128} as a signed byte), so unsigned saturating arithmetic and ordering become
	 * the plain signed saturating ones and the profile holds the raw signed matrix scores. This saves the bias
	 * subtraction per cell and avoids the unsigned operators (whose byte reductions are broken in JDK 25 and crash C2 in
	 * JDK 27). Cells only differ from the C values once a score reaches {@code 255 - bias}, which is reported as
	 * overflow by both.
	 */
	private void swByte(byte[] ref, boolean reverse, int refLen, int readLen, int gapO, int gapE, byte[] profile,
			int terminate, int bias, int maskLen) {
		final int lanes = BS.length();
		final int segLen = (readLen + lanes - 1) / lanes;
		final int stride = segLen * lanes;
		if (bH1.length < stride) {
			int size = Math.max(stride, bH1.length * 2);
			bH1 = new byte[size];
			bH2 = new byte[size];
			bE = new byte[size];
			bHmax = new byte[size];
		}
		Arrays.fill(bH1, 0, stride, ZERO);
		Arrays.fill(bH2, 0, stride, ZERO);
		Arrays.fill(bE, 0, stride, ZERO);
		Arrays.fill(bHmax, 0, stride, ZERO);
		final int[] maxCol = maxColumns(refLen);
		byte[] hStore = bH1, hLoad = bH2;
		final byte[] pvE = bE, hMax = bHmax;

		final ByteVector vZero = ByteVector.broadcast(BS, ZERO);
		// Loaded rather than broadcast: C2 otherwise re-broadcasts them for every use inside the loops.
		Arrays.fill(bGaps, 0, lanes, (byte) gapO);
		Arrays.fill(bGaps, lanes, 2 * lanes, (byte) gapE);
		final ByteVector vGapO = ByteVector.fromArray(BS, bGaps, 0);
		final ByteVector vGapE = ByteVector.fromArray(BS, bGaps, lanes);
		ByteVector vMaxScore = vZero; // highest score of the whole matrix
		ByteVector vMaxMark = vZero; // highest score up to the previous column
		int max = 0;
		int endRef = -1;

		final int begin = reverse ? refLen - 1 : 0, end = reverse ? -1 : refLen, step = reverse ? -1 : 1;
		for (int i = begin; i != end; i += step) {
			ByteVector vF = vZero, vMaxColumn = vZero;
			ByteVector vH = Lanes.shiftBytes(ByteVector.fromArray(BS, hStore, stride - lanes));
			final int pOff = ref[i] * stride;
			byte[] swap = hLoad;
			hLoad = hStore;
			hStore = swap;

			for (int off = 0; off < stride; off += lanes) {
				vH = vH.lanewise(SADD, ByteVector.fromArray(BS, profile, pOff + off));
				ByteVector e = ByteVector.fromArray(BS, pvE, off);
				vH = vH.lanewise(MAX, e).lanewise(MAX, vF);
				vMaxColumn = vMaxColumn.lanewise(MAX, vH);
				vH.intoArray(hStore, off);

				vH = vH.lanewise(SSUB, vGapO);
				e = e.lanewise(SSUB, vGapE).lanewise(MAX, vH);
				e.intoArray(pvE, off);
				vF = vF.lanewise(SSUB, vGapE).lanewise(MAX, vH);

				vH = ByteVector.fromArray(BS, hLoad, off);
			}

			// Lazy-F loop: disallows adjacent insertion and then deletion, so E is not updated (as in SWPS3).
			lazyF: for (int k = 0; k < lanes; ++k) {
				vF = Lanes.shiftBytes(vF);
				for (int off = 0; off < stride; off += lanes) {
					vH = ByteVector.fromArray(BS, hStore, off).lanewise(MAX, vF);
					vMaxColumn = vMaxColumn.lanewise(MAX, vH);
					vH.intoArray(hStore, off);
					vH = vH.lanewise(SSUB, vGapO);
					vF = vF.lanewise(SSUB, vGapE);
					if (!vF.compare(GT, vH).anyTrue()) {
						break lazyF;
					}
				}
			}

			vMaxScore = vMaxScore.lanewise(MAX, vMaxColumn);
			if (!vMaxMark.eq(vMaxScore).allTrue()) {
				vMaxMark = vMaxScore;
				int temp = vMaxScore.reduceLanes(MAX) + 128;
				if (temp > max) {
					max = temp;
					if (max + bias >= 255) {
						break; // overflow
					}
					endRef = i;
					// Keep the column with the best score to locate the end position on the read.
					System.arraycopy(hStore, 0, hMax, 0, stride);
				}
			}
			final int columnMax = vMaxColumn.reduceLanes(MAX) + 128;
			maxCol[i] = columnMax;
			if (columnMax == terminate) {
				break;
			}
		}

		// The smallest query position holding the best score in the stored column.
		int endRead = readLen - 1;
		for (int pos = 0; pos < readLen; pos++) {
			if (hMax[(pos % segLen) * lanes + pos / segLen] + 128 == max) {
				endRead = pos;
				break;
			}
		}
		bestScore = max + bias >= 255 ? 255 : max;
		bestRef = endRef;
		bestRead = endRead;
		secondBest(maxCol, refLen, endRef, maskLen, 1);
	}

	/** 16 bit version of {@link #swByte}, used when the 8 bit kernel saturates. */
	private void swWord(byte[] ref, boolean reverse, int refLen, int readLen, int gapO, int gapE, short[] profile,
			int terminate, int maskLen) {
		final int lanes = WS.length();
		final int segLen = (readLen + lanes - 1) / lanes;
		final int stride = segLen * lanes;
		if (wH1.length < stride) {
			int size = Math.max(stride, wH1.length * 2);
			wH1 = new short[size];
			wH2 = new short[size];
			wE = new short[size];
			wHmax = new short[size];
		} else {
			Arrays.fill(wH1, 0, stride, (short) 0);
			Arrays.fill(wH2, 0, stride, (short) 0);
			Arrays.fill(wE, 0, stride, (short) 0);
			Arrays.fill(wHmax, 0, stride, (short) 0);
		}
		final int[] maxCol = maxColumns(refLen);
		short[] hStore = wH1, hLoad = wH2;
		final short[] pvE = wE, hMax = wHmax;

		final ShortVector vZero = ShortVector.zero(WS);
		Arrays.fill(wGaps, 0, lanes, (short) gapO);
		Arrays.fill(wGaps, lanes, 2 * lanes, (short) gapE);
		final ShortVector vGapO = ShortVector.fromArray(WS, wGaps, 0);
		final ShortVector vGapE = ShortVector.fromArray(WS, wGaps, lanes);
		ShortVector vMaxScore = vZero;
		ShortVector vMaxMark = vZero;
		int max = 0;
		int endRef = 0;

		final int begin = reverse ? refLen - 1 : 0, end = reverse ? -1 : refLen, step = reverse ? -1 : 1;
		for (int i = begin; i != end; i += step) {
			ShortVector vF = vZero, vMaxColumn = vZero;
			ShortVector vH = Lanes.shiftShorts(ShortVector.fromArray(WS, hStore, stride - lanes));
			final int pOff = ref[i] * stride;
			short[] swap = hLoad;
			hLoad = hStore;
			hStore = swap;

			for (int off = 0; off < stride; off += lanes) {
				vH = vH.lanewise(SADD, ShortVector.fromArray(WS, profile, pOff + off));
				ShortVector e = ShortVector.fromArray(WS, pvE, off);
				vH = vH.lanewise(MAX, e).lanewise(MAX, vF);
				vMaxColumn = vMaxColumn.lanewise(MAX, vH);
				vH.intoArray(hStore, off);

				vH = vH.lanewise(SUSUB, vGapO);
				e = e.lanewise(SUSUB, vGapE).lanewise(MAX, vH);
				e.intoArray(pvE, off);
				vF = vF.lanewise(SUSUB, vGapE).lanewise(MAX, vH);

				vH = ShortVector.fromArray(WS, hLoad, off);
			}

			lazyF: for (int k = 0; k < lanes; ++k) {
				vF = Lanes.shiftShorts(vF);
				for (int off = 0; off < stride; off += lanes) {
					vH = ShortVector.fromArray(WS, hStore, off).lanewise(MAX, vF);
					vMaxColumn = vMaxColumn.lanewise(MAX, vH);
					vH.intoArray(hStore, off);
					vH = vH.lanewise(SUSUB, vGapO);
					vF = vF.lanewise(SUSUB, vGapE);
					if (!vF.compare(GT, vH).anyTrue()) {
						break lazyF;
					}
				}
			}

			vMaxScore = vMaxScore.lanewise(MAX, vMaxColumn);
			if (!vMaxMark.eq(vMaxScore).allTrue()) {
				vMaxMark = vMaxScore;
				int temp = vMaxScore.reduceLanes(MAX) & 0xffff;
				if (temp > max) {
					max = temp;
					endRef = i;
					System.arraycopy(hStore, 0, hMax, 0, stride);
				}
			}
			final int columnMax = vMaxColumn.reduceLanes(MAX) & 0xffff;
			maxCol[i] = columnMax;
			if (columnMax == terminate) {
				break;
			}
		}

		int endRead = readLen - 1;
		for (int pos = 0; pos < readLen; pos++) {
			if ((hMax[(pos % segLen) * lanes + pos / segLen] & 0xffff) == max) {
				endRead = pos;
				break;
			}
		}
		bestScore = max;
		bestRef = endRef;
		bestRead = endRead;
		secondBest(maxCol, refLen, endRef, maskLen, 0);
	}

	private int[] maxColumns(int refLen) {
		if (maxColumn.length < refLen) {
			maxColumn = new int[Math.max(refLen, maxColumn.length * 2)];
		} else {
			Arrays.fill(maxColumn, 0, refLen, 0);
		}
		return maxColumn;
	}

	/**
	 * Highest column maximum outside {@code [endRef - maskLen, endRef + maskLen]}. The C byte and word kernels differ
	 * by one in where the right hand side starts, kept here through {@code rightOffset}.
	 */
	private void secondBest(int[] maxCol, int refLen, int endRef, int maskLen, int rightOffset) {
		int score = 0, pos = 0;
		int edge = Math.max(endRef - maskLen, 0);
		for (int i = 0; i < edge; i++) {
			if (maxCol[i] > score) {
				score = maxCol[i];
				pos = i;
			}
		}
		edge = Math.min(endRef + maskLen, refLen);
		for (int i = edge + rightOffset; i < refLen; i++) {
			if (maxCol[i] > score) {
				score = maxCol[i];
				pos = i;
			}
		}
		secondScore = score;
		secondRef = pos;
	}

	/**
	 * Banded global alignment of the already located best local alignment region, with traceback. Returns the packed
	 * cigar, or {@code null} when the traceback fails.
	 */
	private int[] bandedSw(byte[] ref, int refOff, byte[] read, int readOff, int refLen, int readLen, int score,
			int gapO, int gapE, int bandWidth, byte[] matT, int n) {
		int max = 0, maxI = 0, maxJ = 0;
		final int len = Math.max(refLen, readLen);
		int widthD;
		do {
			final int width = bandWidth * 2 + 3;
			widthD = bandWidth * 2 + 1;
			if (width >= hB.length) {
				int size = Integer.highestOneBit(width) << 1;
				hB = Arrays.copyOf(hB, size);
				eB = Arrays.copyOf(eB, size);
				hC = Arrays.copyOf(hC, size);
			}
			final long dirSize = (long) widthD * readLen;
			if (dirSize >= direction.length) {
				if (dirSize >= Integer.MAX_VALUE - 8) {
					throw new IllegalStateException("alignment band too large: " + dirSize);
				}
				direction = new byte[(int) Math.min(Integer.MAX_VALUE - 8, Long.highestOneBit(dirSize) << 1)];
			}
			final int[] hb = hB, eb = eB, hc = hC;
			final byte[] dir = direction;
			Arrays.fill(hb, 1, width - 1, 0);
			for (int i = 0; i < readLen; i++) {
				final int end = Math.min(refLen - 1, i + bandWidth);
				final int edge = Math.min(end + 1, width - 1);
				final int x = Math.max(i - bandWidth, 0); // band start of this row (= first j)
				final int delta = x - Math.max(i - 1 - bandWidth, 0); // shift against the previous row's band
				final int uEnd = end - x + 1;
				final int rowOff = read[readOff + i] * n;
				hb[0] = hb[edge] = hc[0] = 0;
				eb[0] = eb[edge] = NEG_INF;
				final long best = bandRow(hb, eb, hc, dir, widthD * i - 1, ref, refOff + x - 1, matT, rowOff, uEnd,
						delta, i == 0, gapO, gapE);
				if ((int) (best >>> 32) > max) {
					max = (int) (best >>> 32);
					maxI = i;
					maxJ = x + (int) best - 1;
				}
				if (uEnd > 0) {
					System.arraycopy(hc, 1, hb, 1, uEnd);
				}
			}
			bandWidth *= 2;
		} while (max < score && bandWidth <= len);
		bandWidth /= 2;
		widthD = bandWidth * 2 + 1;

		// Trace back.
		int i = maxI, j = maxJ;
		int count = 0, l = 0;
		int op = Cigar.M, prevOp = Cigar.M;
		int state = 2; // 0: e, 1: f, 2: h
		int line = widthD * i;
		int[] c = cigarBuf;
		final byte[] dir = direction;
		while (i >= 0 && j > 0) {
			final int x = Math.max(i - bandWidth, 0);
			final int d = dir[line + j - x];
			final int code = switch (state) {
				case 0 -> (d & DIR_E_OPEN) != 0 ? 3 : 2;
				case 1 -> (d & DIR_F_OPEN) != 0 ? 5 : 4;
				default -> d & 7;
			};
			switch (code) {
				case 1 -> {
					--i;
					--j;
					state = 2;
					line -= widthD;
					op = Cigar.M;
				}
				case 2 -> {
					--i;
					state = 0;
					line -= widthD;
					op = Cigar.I;
				}
				case 3 -> {
					--i;
					state = 2;
					line -= widthD;
					op = Cigar.I;
				}
				case 4 -> {
					--j;
					state = 1;
					op = Cigar.D;
				}
				case 5 -> {
					--j;
					state = 2;
					op = Cigar.D;
				}
				default -> {
					return null;
				}
			}
			if (op == prevOp) {
				++count;
			} else {
				if (l + 2 >= c.length) {
					c = cigarBuf = Arrays.copyOf(c, c.length * 2);
				}
				c[l++] = Cigar.toInt(count, prevOp);
				prevOp = op;
				count = 1;
			}
		}
		if (l + 2 >= c.length) {
			c = cigarBuf = Arrays.copyOf(c, c.length * 2);
		}
		if (op == Cigar.M) {
			c[l++] = Cigar.toInt(count + 1, op);
		} else {
			c[l++] = Cigar.toInt(count, op);
			c[l++] = Cigar.toInt(1, Cigar.M);
		}
		int[] result = new int[l];
		for (int k = 0; k < l; k++) {
			result[k] = c[l - 1 - k];
		}
		return result;
	}

	/**
	 * One row of the banded alignment. Returns the row's highest H in the upper 32 bits and the first band index
	 * reaching it in the lower 32 bits. Kept separate from {@link #bandedSw} so the JIT has registers to spare.
	 *
	 * @param matT    transposed substitution matrix, {@code matT[rowOff + refLetter]} scores this query letter
	 * @param delta   1 when this row's band starts one column further than the previous row's, else 0
	 */
	private static long bandRow(int[] hb, int[] eb, int[] hc, byte[] dir, int lineBase, byte[] ref, int refBase,
			byte[] matT, int rowOff, int uEnd, int delta, boolean first, int gapO, int gapE) {
		int f = NEG_INF;
		int rowMax = Integer.MIN_VALUE, rowMaxU = 0;
		// u is the band index of column j = x + u - 1, e the band index of the same column in the previous row.
		// The direction choices are random per cell, so they are computed from sign bits instead of branches (C2
		// keeps them as mispredicting jumps). No difference overflows: all values are >= NEG_INF - 255.
		for (int u = 1; u <= uEnd; u++) {
			final int e = u + delta;
			int temp1 = first ? -gapO : hb[e] - gapO;
			int temp2 = first ? NEG_INF : eb[e] - gapE;
			final int eOpen = (temp2 - temp1) >>> 31; // temp1 > temp2
			final int eu = Math.max(temp1, temp2);
			eb[u] = eu;

			temp1 = hc[u - 1] - gapO;
			temp2 = f - gapE;
			final int fOpen = (temp2 - temp1) >>> 31;
			f = Math.max(temp1, temp2);

			final int e1 = Math.max(eu, 0);
			final int f1 = Math.max(f, 0);
			temp1 = Math.max(e1, f1);
			temp2 = hb[e - 1] + matT[rowOff + ref[refBase + u]];
			final int h = Math.max(temp1, temp2);
			hc[u] = h;
			if (h > rowMax) {
				rowMax = h;
				rowMaxU = u;
			}
			// H from E (C codes 2/3) when e1 > f1, else from F (4/5); diagonal (1) when temp1 <= temp2.
			final int fromE = (f1 - e1) >>> 31;
			final int gap = 4 - 2 * fromE + fOpen + fromE * (eOpen - fOpen);
			final int diagonal = ((temp2 - temp1) >>> 31) ^ 1;
			final int dirH = gap + diagonal * (1 - gap);
			dir[lineBase + u] = (byte) (dirH | eOpen << 3 | fOpen << 4);
		}
		return ((long) rowMax << 32) | rowMaxU;
	}

	/**
	 * Re-scores a banded path. A path from a too narrow band can run past the aligned region (the C code then reads
	 * out of bounds); that is reported as {@link Integer#MIN_VALUE} so the caller retries with the full band.
	 */
	private static int cigarScore(int[] path, byte[] ref, int refPos, int refLen, byte[] read, int readPos,
			int readLen, byte[] mat, int n, int gapO, int gapE) {
		final int refEnd = refPos + refLen, readEnd = readPos + readLen;
		int score = 0;
		for (int c : path) {
			final int len = Cigar.length(c);
			final int op = Cigar.op(c);
			if (op == Cigar.M) {
				if (refPos + len > refEnd || readPos + len > readEnd) {
					return Integer.MIN_VALUE;
				}
				for (int j = 0; j < len; ++j) {
					score += mat[ref[refPos++] * n + read[readPos++]];
				}
			} else {
				score -= gapO + (len > 1 ? (len - 1) * gapE : 0);
				if (op == Cigar.I) {
					readPos += len;
				} else if (op == Cigar.D) {
					refPos += len;
				}
			}
		}
		return score;
	}
}
