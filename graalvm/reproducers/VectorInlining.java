import static jdk.incubator.vector.VectorOperators.AND;
import static jdk.incubator.vector.VectorOperators.GT;
import static jdk.incubator.vector.VectorOperators.LSHL;
import static jdk.incubator.vector.VectorOperators.LSHR;
import static jdk.incubator.vector.VectorOperators.MAX;
import static jdk.incubator.vector.VectorOperators.OR;
import static jdk.incubator.vector.VectorOperators.SADD;
import static jdk.incubator.vector.VectorOperators.SSUB;

import java.util.Arrays;
import java.util.Random;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

/**
 * On GraalVM, when a method that takes or returns a Vector API vector is not inlined into the hot loop, every vector
 * operation of that loop runs its Java fallback ({@code ByteVector.bOpTemplate} etc.), about 30x slower than C2 in the
 * same situation and about 140x slower than with the method inlined.
 * <p>
 * The kernel is the 8-bit striped Smith-Waterman column loop of cssw4j (a Java port of the SSW library), reduced to
 * the score. {@link #shiftBytes} shifts the lanes of a vector by one byte; it is called twice per target column.
 * <p>
 * Run it as is, then with {@code shiftBytes} kept out of line:
 *
 * <pre>
 * java --add-modules jdk.incubator.vector VectorInlining.java
 * java -XX:CompileCommand=quiet -XX:CompileCommand=dontinline,VectorInlining::shiftBytes \
 *      --add-modules jdk.incubator.vector VectorInlining.java
 * </pre>
 *
 * In a real program (a database search with cssw4j) nothing forces this: Graal's inliner left the equivalent helper
 * out of line by itself ("not-within-budget" in {@code -Djdk.graal.Log=Inlining}) once the calling loop had grown, and
 * the search became about 70x slower.
 */
public class VectorInlining {
	static final VectorSpecies<Byte> BS = ByteVector.SPECIES_PREFERRED;
	static final VectorSpecies<Integer> INT = IntVector.SPECIES_PREFERRED;
	static final VectorShuffle<Integer> PREVIOUS = VectorShuffle.fromOp(INT, i -> Math.max(i - 1, 0));
	static final int[] NOT_FIRST = new int[INT.length()], BYTE_FILL = new int[INT.length()];
	/** Score 0 in excess-128 form. */
	static final byte ZERO = Byte.MIN_VALUE;

	static {
		Arrays.fill(NOT_FIRST, 1, NOT_FIRST.length, -1);
		BYTE_FILL[0] = 0x80;
	}

	byte[] h1 = new byte[0], h2 = new byte[0], e = new byte[0];
	final byte[] gaps = new byte[2 * BS.length()];

	/** Lane k of the result is lane k - 1 of v; lane 0 becomes ZERO. */
	static ByteVector shiftBytes(ByteVector v) {
		IntVector l = (IntVector) v.reinterpretAsInts();
		IntVector carry = l.rearrange(PREVIOUS).lanewise(LSHR, 24).lanewise(AND, IntVector.fromArray(INT, NOT_FIRST, 0));
		return l.lanewise(LSHL, 8).lanewise(OR, carry).lanewise(OR, IntVector.fromArray(INT, BYTE_FILL, 0))
				.reinterpretAsBytes();
	}

	/** Best local alignment score of the striped profile against target (forward or reversed). */
	int kernel(byte[] target, boolean reverse, int queryLength, int gapO, int gapE, byte[] profile) {
		final int lanes = BS.length();
		final int stride = (queryLength + lanes - 1) / lanes * lanes;
		if (h1.length < stride) {
			h1 = new byte[stride];
			h2 = new byte[stride];
			e = new byte[stride];
		}
		Arrays.fill(h1, 0, stride, ZERO);
		Arrays.fill(h2, 0, stride, ZERO);
		Arrays.fill(e, 0, stride, ZERO);
		byte[] hStore = h1, hLoad = h2;
		final byte[] pvE = e;
		final ByteVector vZero = ByteVector.broadcast(BS, ZERO);
		Arrays.fill(gaps, 0, lanes, (byte) gapO);
		Arrays.fill(gaps, lanes, 2 * lanes, (byte) gapE);
		final ByteVector vGapO = ByteVector.fromArray(BS, gaps, 0);
		final ByteVector vGapE = ByteVector.fromArray(BS, gaps, lanes);
		ByteVector vMaxScore = vZero, vMaxMark = vZero;
		int max = 0;
		final int n = target.length;
		final int begin = reverse ? n - 1 : 0, end = reverse ? -1 : n, step = reverse ? -1 : 1;
		for (int i = begin; i != end; i += step) {
			ByteVector vF = vZero, vMaxColumn = vZero;
			ByteVector vH = shiftBytes(ByteVector.fromArray(BS, hStore, stride - lanes));
			final int pOff = target[i] * stride;
			byte[] swap = hLoad;
			hLoad = hStore;
			hStore = swap;
			for (int off = 0; off < stride; off += lanes) {
				vH = vH.lanewise(SADD, ByteVector.fromArray(BS, profile, pOff + off));
				ByteVector ve = ByteVector.fromArray(BS, pvE, off);
				vH = vH.lanewise(MAX, ve).lanewise(MAX, vF);
				vMaxColumn = vMaxColumn.lanewise(MAX, vH);
				vH.intoArray(hStore, off);
				vH = vH.lanewise(SSUB, vGapO);
				ve = ve.lanewise(SSUB, vGapE).lanewise(MAX, vH);
				ve.intoArray(pvE, off);
				vF = vF.lanewise(SSUB, vGapE).lanewise(MAX, vH);
				vH = ByteVector.fromArray(BS, hLoad, off);
			}
			lazyF: for (int k = 0; k < lanes; ++k) {
				vF = shiftBytes(vF);
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
				max = Math.max(max, vMaxScore.reduceLanes(MAX) + 128);
				if (max >= 250) {
					break; // overflow, not reached in this test
				}
			}
		}
		return max;
	}

	/** The hot loop over all targets, as in a database search. */
	static long scan(VectorInlining a, byte[][] targets, int queryLength, byte[] profile) {
		long sum = 0;
		for (byte[] t : targets) {
			sum += a.kernel(t, false, queryLength, 12, 1, profile);
		}
		return sum;
	}

	public static void main(String[] args) {
		Random rnd = new Random(1);
		int queryLength = 150, stride = (queryLength + BS.length() - 1) / BS.length() * BS.length();
		byte[] query = new byte[queryLength];
		for (int i = 0; i < queryLength; i++) {
			query[i] = (byte) rnd.nextInt(20);
		}
		// Striped profile (SSW layout): element (segment s, lane l) scores query position l * segments + s;
		// +5 for a match, -1..-3 otherwise and for the padding.
		byte[] profile = new byte[20 * stride];
		int lanes = BS.length(), segments = stride / lanes;
		for (int letter = 0; letter < 20; letter++) {
			for (int s = 0; s < segments; s++) {
				for (int l = 0; l < lanes; l++) {
					int pos = l * segments + s;
					boolean match = pos < queryLength && query[pos] == letter;
					profile[letter * stride + s * lanes + l] = (byte) (match ? 5 : -1 - rnd.nextInt(3));
				}
			}
		}
		// 1 in 20 targets contains a 25-residue piece of the query: a hit.
		byte[][] targets = new byte[20000][];
		for (int i = 0; i < targets.length; i++) {
			targets[i] = new byte[50 + rnd.nextInt(600)];
			for (int j = 0; j < targets[i].length; j++) {
				targets[i][j] = (byte) rnd.nextInt(20);
			}
			if (i % 20 == 0) {
				System.arraycopy(query, 40, targets[i], 10, 25);
			}
		}
		VectorInlining a = new VectorInlining();
		System.out.println(System.getProperty("java.vm.name") + " " + System.getProperty("java.vm.version") + ", " + BS);
		for (int round = 0; round < 8; round++) {
			long t0 = System.nanoTime();
			long sum = scan(a, targets, queryLength, profile);
			System.out.printf("round %2d: %7.0f ms (checksum %d)%n", round, (System.nanoTime() - t0) / 1e6, sum);
		}
	}
}
