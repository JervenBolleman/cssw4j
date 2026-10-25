package swiss.sib.cssw4j.tornado;

import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.annotations.Parallel;

/**
 * Smith-Waterman scoring kernel for TornadoVM.
 * <p>
 * The striped SIMD layout of the CPU code does not map to GPU lanes: it packs one alignment into the lanes of one
 * vector register. A GPU has thousands of threads instead, so this kernel uses inter-sequence parallelism, as
 * CUDASW++ does: every work-item aligns one query against one target with a plain scalar affine-gap
 * (Gotoh) recurrence. The kernel is ordinary Java, so it also runs sequentially on the JVM, which the tests use to
 * check the algorithm without a device.
 * <p>
 * For every pair (query q, target t), with {@code p = q * numTargets + t}, it computes the optimal local alignment
 * score and its end positions on the target and the query, with the tie-breaking of the SSW C library: the first
 * target column that reaches the maximum, and in that column the smallest query position (0 and 0 when the score is
 * 0). Begin positions and the CIGAR are not computed here; run the CPU aligner on the hits that matter.
 * <p>
 * Scores are clamped at 0 for the gap states as well, like the saturating arithmetic of the SSW kernels. As in SSW,
 * the scores equal the exact Smith-Waterman scores for gapOpen &gt; gapExtend.
 */
public final class SwKernel {

	private SwKernel() {
	}

	/**
	 * @param params     {@code [gapOpen, gapExtend, matrixSize, jStep, pStep]}; the scratch index of query position
	 *                   {@code j} of pair {@code p} is {@code j * jStep + p * pStep} (see {@link Layout})
	 * @param matrix     row major substitution matrix, row = target letter, column = query letter
	 * @param queries    all queries, encoded as matrix indices, one after another
	 * @param qOff       start of every query in {@code queries}, plus the total length at the end
	 * @param targets    all targets, encoded, one after another
	 * @param tOff       start of every target in {@code targets}, plus the total length at the end
	 * @param h          scratch, {@code (maxQueryLength + 1) * numPairs} ints
	 * @param e          scratch, same size as {@code h}
	 * @param scores     out: best score of every pair
	 * @param refEnds    out: 0-based end position on the target of the best alignment
	 * @param readEnds   out: 0-based end position on the query of the best alignment
	 */
	public static void scorePairs(int[] params, int[] matrix, byte[] queries, int[] qOff, byte[] targets, int[] tOff,
			int[] h, int[] e, int[] scores, int[] refEnds, int[] readEnds) {
		final int nPairs = scores.length;
		final int nTargets = tOff.length - 1;
		final int gapO = params[0];
		final int gapE = params[1];
		final int n = params[2];
		final int jStep = params[3];
		final int pStep = params[4];
		for (@Parallel int p = 0; p < nPairs; p++) {
			final int q = p / nTargets;
			final int t = p - q * nTargets;
			final int qStart = qOff[q];
			final int qLen = qOff[q + 1] - qStart;
			final int tStart = tOff[t];
			final int tEnd = tOff[t + 1];
			for (int j = 0; j <= qLen; j++) {
				h[j * jStep + p * pStep] = 0;
				e[j * jStep + p * pStep] = 0;
			}
			int best = 0;
			int bestRef = 0; // like the C library, the end positions of a pair without any match are 0
			int bestRead = 0;
			for (int i = tStart; i < tEnd; i++) {
				final int row = targets[i] * n;
				int diag = 0; // H of the previous column, one row up
				int up = 0; // H of this column, one row up
				int f = 0; // vertical gap state of this column
				for (int j = 1; j <= qLen; j++) {
					final int idx = j * jStep + p * pStep;
					final int left = h[idx]; // H of the previous column, this row
					final int gapH = Math.max(Math.max(e[idx] - gapE, left - gapO), 0);
					f = Math.max(Math.max(f - gapE, up - gapO), 0);
					final int score = matrix[row + queries[qStart + j - 1]];
					final int cell = Math.max(Math.max(diag + score, 0), Math.max(gapH, f));
					h[idx] = cell;
					e[idx] = gapH;
					diag = left;
					up = cell;
					if (cell > best) {
						best = cell;
						bestRef = i - tStart;
						bestRead = j - 1;
					}
				}
			}
			scores[p] = best;
			refEnds[p] = bestRef;
			readEnds[p] = bestRead;
		}
	}

	/**
	 * Same result as {@link #scorePairs}, with less global memory traffic: the query is processed in strips of four
	 * rows whose H and E values stay in registers while the strip walks the whole target. Only the bottom row of a strip
	 * (its H and F values, one per target position) goes through global memory, to the next strip. That is 16 bytes per
	 * four cells instead of 16 bytes per cell.
	 * <p>
	 * The scratch arrays hold one value per target position (not query position): {@code bH} and {@code bF} must have
	 * {@code (maxTargetLength + 1) * numPairs} ints, indexed {@code i * params[3] + p * params[4]}.
	 * <p>
	 * The strips change the order in which cells are visited, so ties are resolved explicitly: among cells with the best
	 * score the one with the smallest target position wins, then the smallest query position, as in {@link #scorePairs}.
	 * <p>
	 * The query letters are passed as ints ({@link Sequences#residuesAsInts()}): with byte letters TornadoVM 7.3.1
	 * generated {@code letter + row} with the byte type of the letter, which truncated the matrix index on every
	 * backend.
	 */
	public static void scorePairsStrip4(int[] params, int[] matrix, int[] queries, int[] qOff, byte[] targets,
			int[] tOff, int[] bH, int[] bF, int[] scores, int[] refEnds, int[] readEnds) {
		for (@Parallel int p = 0; p < scores.length; p++) {
			strip4Pair(p, params, matrix, queries, qOff, targets, tOff, bH, bF, scores, refEnds, readEnds);
		}
	}

	/** Largest substitution matrix (n * n entries) supported by {@link #scorePairsStrip4Local}. */
	public static final int MAX_MATRIX_CELLS = 1024;

	/**
	 * {@link #scorePairsStrip4} with the substitution matrix copied into the work-group's local (shared) memory first.
	 * Every cell reads the matrix at a position that depends on the target letter, which differs between neighbouring
	 * work-items; from global memory those scattered reads take more than half of the kernel time on the Quadro M1000M.
	 * Needs a launch with a local work size (see {@link TornadoScorer}), and a matrix of at most
	 * {@link #MAX_MATRIX_CELLS} entries.
	 */
	public static void scorePairsStrip4Local(KernelContext context, int[] params, int[] matrix, int[] queries,
			int[] qOff, byte[] targets, int[] tOff, int[] bH, int[] bF, int[] scores, int[] refEnds, int[] readEnds) {
		final int[] local = context.allocateIntLocalArray(MAX_MATRIX_CELLS);
		final int cells = params[2] * params[2];
		for (int k = context.localIdx; k < cells; k += context.localGroupSizeX) {
			local[k] = matrix[k];
		}
		context.localBarrier();
		final int p = context.globalIdx;
		if (p < scores.length) {
			strip4Pair(p, params, local, queries, qOff, targets, tOff, bH, bF, scores, refEnds, readEnds);
		}
	}

	/** The strip-of-four recurrence for pair {@code p}; shared by the strip kernels. */
	private static void strip4Pair(int p, int[] params, int[] matrix, int[] queries, int[] qOff, byte[] targets,
			int[] tOff, int[] bH, int[] bF, int[] scores, int[] refEnds, int[] readEnds) {
		final int nTargets = tOff.length - 1;
		final int gapO = params[0];
		final int gapE = params[1];
		final int n = params[2];
		final int iStep = params[3];
		final int pStep = params[4];
		final int q = p / nTargets;
		final int t = p - q * nTargets;
		final int qStart = qOff[q];
		final int qLen = qOff[q + 1] - qStart;
		final int tStart = tOff[t];
		final int tLen = tOff[t + 1] - tStart;
		for (int i = 0; i <= tLen; i++) {
			bH[i * iStep + p * pStep] = 0; // row 0: H and F are 0
			bF[i * iStep + p * pStep] = 0;
		}
		int best = 0;
		int bestRef = 0;
		int bestRead = 0;
		for (int j0 = 0; j0 < qLen; j0 += 4) {
			// Query letters of the strip. Rows past the end of the query are padding and never reported; they reuse
			// the last letter so that every load stays inside the array (a conditional load may be compiled into
			// an unconditional one on the device).
			final int last = qStart + qLen - 1;
			final int q0 = queries[qStart + j0];
			final int q1 = queries[Math.min(qStart + j0 + 1, last)];
			final int q2 = queries[Math.min(qStart + j0 + 2, last)];
			final int q3 = queries[Math.min(qStart + j0 + 3, last)];
			int h0 = 0, h1 = 0, h2 = 0, h3 = 0; // H of the previous column (column 0 is 0)
			int e0 = 0, e1 = 0, e2 = 0, e3 = 0; // E of the previous column
			int topPrev = 0; // H of the previous column in the row above the strip
			for (int i = 1; i <= tLen; i++) {
				final int bidx = i * iStep + p * pStep;
				final int top = bH[bidx]; // H of this column in the row above the strip
				int f = bF[bidx]; // F of this column in the row above the strip
				final int row = targets[tStart + i - 1] * n;

				int g = Math.max(Math.max(e0 - gapE, h0 - gapO), 0);
				f = Math.max(Math.max(f - gapE, top - gapO), 0);
				final int c0 = Math.max(Math.max(topPrev + matrix[row + q0], 0), Math.max(g, f));
				final int d1 = h0;
				h0 = c0;
				e0 = g;

				g = Math.max(Math.max(e1 - gapE, h1 - gapO), 0);
				f = Math.max(Math.max(f - gapE, c0 - gapO), 0);
				final int c1 = Math.max(Math.max(d1 + matrix[row + q1], 0), Math.max(g, f));
				final int d2 = h1;
				h1 = c1;
				e1 = g;

				g = Math.max(Math.max(e2 - gapE, h2 - gapO), 0);
				f = Math.max(Math.max(f - gapE, c1 - gapO), 0);
				final int c2 = Math.max(Math.max(d2 + matrix[row + q2], 0), Math.max(g, f));
				final int d3 = h2;
				h2 = c2;
				e2 = g;

				g = Math.max(Math.max(e3 - gapE, h3 - gapO), 0);
				f = Math.max(Math.max(f - gapE, c2 - gapO), 0);
				final int c3 = Math.max(Math.max(d3 + matrix[row + q3], 0), Math.max(g, f));
				h3 = c3;
				e3 = g;

				bH[bidx] = c3; // the bottom row of the strip, for the next strip
				bF[bidx] = f;
				topPrev = top;

				// Best cell, lexicographic tie-break on (target position, query position); padding rows excluded.
				final int ref = i - 1;
				if (c0 > best || c0 == best && c0 > 0 && (ref < bestRef || ref == bestRef && j0 < bestRead)) {
					best = c0;
					bestRef = ref;
					bestRead = j0;
				}
				if (j0 + 1 < qLen && (c1 > best
						|| c1 == best && c1 > 0 && (ref < bestRef || ref == bestRef && j0 + 1 < bestRead))) {
					best = c1;
					bestRef = ref;
					bestRead = j0 + 1;
				}
				if (j0 + 2 < qLen && (c2 > best
						|| c2 == best && c2 > 0 && (ref < bestRef || ref == bestRef && j0 + 2 < bestRead))) {
					best = c2;
					bestRef = ref;
					bestRead = j0 + 2;
				}
				if (j0 + 3 < qLen && (c3 > best
						|| c3 == best && c3 > 0 && (ref < bestRef || ref == bestRef && j0 + 3 < bestRead))) {
					best = c3;
					bestRef = ref;
					bestRead = j0 + 3;
				}
			}
		}
		scores[p] = best;
		refEnds[p] = bestRef;
		readEnds[p] = bestRead;
	}
}
