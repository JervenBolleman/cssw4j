package swiss.sib.cssw4j;

import java.util.Arrays;

/**
 * Striped query profile, the Java counterpart of {@code ssw_init}.
 * <p>
 * Build it once per query and align it against many references. Instances are immutable and may be shared between
 * threads (each thread using its own {@link Aligner}).
 */
public final class Profile {
	/** Which kernels the profile is prepared for, mirroring the {@code score_size} argument of {@code ssw_init}. */
	public enum ScoreSize {
		/** 8 bit kernel only; use when the best score is surely below 255. */
		BYTE,
		/** 16 bit kernel only; use when the best score is surely at least 255. */
		WORD,
		/** 8 bit kernel, falling back to 16 bit when it overflows. */
		BOTH
	}

	final byte[] read;
	final int readLen;
	final ScoringMatrix matrix;
	final int bias;
	final byte[] profileByte;
	final short[] profileWord;

	private Profile(byte[] read, int readLen, ScoringMatrix matrix, ScoreSize size) {
		if (readLen <= 0 || readLen > read.length) {
			throw new IllegalArgumentException("invalid query length " + readLen);
		}
		this.read = read;
		this.readLen = readLen;
		this.matrix = matrix;
		if (size != ScoreSize.WORD) {
			this.bias = -matrix.minScore();
			this.profileByte = buildByte(read, readLen, matrix.matrix(), matrix.size(), null);
		} else {
			this.bias = 0;
			this.profileByte = null;
		}
		this.profileWord = size != ScoreSize.BYTE ? buildWord(read, readLen, matrix.matrix(), matrix.size(), null) : null;
	}

	/**
	 * Profile for a query already encoded with {@link ScoringMatrix#encode}. The array is retained, not copied.
	 */
	public static Profile of(byte[] read, ScoringMatrix matrix) {
		return new Profile(read, read.length, matrix, ScoreSize.BOTH);
	}

	/**
	 * Profile for the first {@code readLen} letters of an encoded query. The array is retained, not copied.
	 */
	public static Profile of(byte[] read, int readLen, ScoringMatrix matrix, ScoreSize size) {
		return new Profile(read, readLen, matrix, size);
	}

	/** Encodes and profiles a query given as letters. */
	public static Profile of(CharSequence read, ScoringMatrix matrix) {
		return of(matrix.encode(read), matrix);
	}

	public int length() {
		return readLen;
	}

	public ScoringMatrix matrix() {
		return matrix;
	}

	/**
	 * Striped 8 bit profile: for every letter {@code nt}, {@code segLen} vectors where lane {@code k} of segment
	 * {@code i} holds the score of query position {@code i + k * segLen} (0 for the padding past the query end). Unlike
	 * the C code the scores are not biased, the 8 bit kernel works in excess-128 form instead.
	 */
	static byte[] buildByte(byte[] read, int readLen, byte[] mat, int n, byte[] reuse) {
		final int lanes = Lanes.BYTE.length();
		final int stride = (readLen + lanes - 1) / lanes * lanes;
		final int size = n * stride;
		byte[] t = reuse != null && reuse.length >= size ? reuse : new byte[size];
		final byte[] striped = stripe(read, readLen, lanes, n);
		final byte[] row = new byte[n + 1]; // row[n] = 0 scores the padding
		for (int nt = 0, k = 0; nt < n; nt++) {
			System.arraycopy(mat, nt * n, row, 0, n);
			for (int p = 0; p < stride; p++) {
				t[k++] = row[striped[p]];
			}
		}
		return t;
	}

	/** Striped 16 bit profile, same layout as {@link #buildByte}. */
	static short[] buildWord(byte[] read, int readLen, byte[] mat, int n, short[] reuse) {
		final int lanes = Lanes.SHORT.length();
		final int stride = (readLen + lanes - 1) / lanes * lanes;
		final int size = n * stride;
		short[] t = reuse != null && reuse.length >= size ? reuse : new short[size];
		final byte[] striped = stripe(read, readLen, lanes, n);
		final short[] row = new short[n + 1];
		for (int nt = 0, k = 0; nt < n; nt++) {
			for (int c = 0; c < n; c++) {
				row[c] = mat[nt * n + c];
			}
			for (int p = 0; p < stride; p++) {
				t[k++] = row[striped[p]];
			}
		}
		return t;
	}

	/**
	 * The query letters in striped order: position {@code s * lanes + k} holds letter {@code k * segLen + s}, or
	 * {@code pad} past the end of the query.
	 */
	private static byte[] stripe(byte[] read, int readLen, int lanes, int pad) {
		final int segLen = (readLen + lanes - 1) / lanes;
		final byte[] striped = new byte[segLen * lanes];
		Arrays.fill(striped, (byte) pad);
		for (int j = 0, k = 0, seg = 0; j < readLen; j++) {
			striped[seg * lanes + k] = read[j];
			if (++seg == segLen) {
				seg = 0;
				k++;
			}
		}
		return striped;
	}
}
