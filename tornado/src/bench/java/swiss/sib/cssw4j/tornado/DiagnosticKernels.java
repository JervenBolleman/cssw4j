package swiss.sib.cssw4j.tornado;

import uk.ac.manchester.tornado.api.annotations.Parallel;

/**
 * Diagnostic kernels for the GPU investigation. They do not compute alignment scores and must not be used for results.
 */
public final class DiagnosticKernels {
	private DiagnosticKernels() {
	}

	/**
	 * {@link SwKernel#scorePairsStrip4} with the substitution-matrix lookup replaced by a letter comparison (+5 for
	 * equal letters, -2 otherwise). Same recurrence, same loads of H, F and the letters, but no read of the matrix: the
	 * difference in run time is the cost of the matrix lookups. The {@code matrix} argument is unused and kept for an
	 * identical task-graph shape.
	 */
	public static void strip4LetterScore(int[] params, int[] matrix, int[] queries, int[] qOff, byte[] targets,
			int[] tOff, int[] bH, int[] bF, int[] scores, int[] refEnds, int[] readEnds) {
		final int nPairs = scores.length;
		final int nTargets = tOff.length - 1;
		final int gapO = params[0];
		final int gapE = params[1];
		final int iStep = params[3];
		final int pStep = params[4];
		for (@Parallel int p = 0; p < nPairs; p++) {
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
					final int tl = targets[tStart + i - 1];

					int g = Math.max(Math.max(e0 - gapE, h0 - gapO), 0);
					f = Math.max(Math.max(f - gapE, top - gapO), 0);
					final int c0 = Math.max(Math.max(topPrev + (tl == q0 ? 5 : -2), 0), Math.max(g, f));
					final int d1 = h0;
					h0 = c0;
					e0 = g;

					g = Math.max(Math.max(e1 - gapE, h1 - gapO), 0);
					f = Math.max(Math.max(f - gapE, c0 - gapO), 0);
					final int c1 = Math.max(Math.max(d1 + (tl == q1 ? 5 : -2), 0), Math.max(g, f));
					final int d2 = h1;
					h1 = c1;
					e1 = g;

					g = Math.max(Math.max(e2 - gapE, h2 - gapO), 0);
					f = Math.max(Math.max(f - gapE, c1 - gapO), 0);
					final int c2 = Math.max(Math.max(d2 + (tl == q2 ? 5 : -2), 0), Math.max(g, f));
					final int d3 = h2;
					h2 = c2;
					e2 = g;

					g = Math.max(Math.max(e3 - gapE, h3 - gapO), 0);
					f = Math.max(Math.max(f - gapE, c2 - gapO), 0);
					final int c3 = Math.max(Math.max(d3 + (tl == q3 ? 5 : -2), 0), Math.max(g, f));
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
}
