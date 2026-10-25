package swiss.sib.cssw4j.tornado;

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
}
