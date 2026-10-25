package swiss.sib.cssw4j.tornado;

import uk.ac.manchester.tornado.api.annotations.Parallel;

/** A kernel that only writes the three result arrays, to measure everything around the alignment work. */
public final class OverheadKernel {
	private OverheadKernel() {
	}

	/** Has the same arguments as {@link SwKernel#scorePairs}, so the task graph transfers the same data. */
	public static void clear(int[] params, int[] matrix, byte[] queries, int[] qOff, byte[] targets, int[] tOff, int[] h,
			int[] e, int[] scores, int[] refEnds, int[] readEnds) {
		for (@Parallel int p = 0; p < scores.length; p++) {
			scores[p] = 0;
			refEnds[p] = 0;
			readEnds[p] = 0;
		}
	}
}
