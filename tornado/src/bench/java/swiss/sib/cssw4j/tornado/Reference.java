package swiss.sib.cssw4j.tornado;

import java.util.Arrays;

/**
 * Checks benchmark results against {@link SwKernel#scorePairs} run as plain Java on the CPU, which matches the C SSW
 * library on the test dataset (see Check). Benchmarks call it once in their setup, so a faster but wrong kernel fails the
 * run instead of producing a number.
 */
final class Reference {
	private Reference() {
	}

	static void check(Matrix matrix, Sequences queries, Sequences targets, TornadoScorer scorer) {
		int pairs = queries.count() * targets.count();
		int[] steps = Layout.INTERLEAVED.steps(pairs, queries.maxLength);
		int[] params = { Workloads.GAP_OPEN, Workloads.GAP_EXTEND, matrix.size(), steps[0], steps[1] };
		int[] h = new int[(queries.maxLength + 1) * pairs], e = new int[h.length];
		int[] scores = new int[pairs], refEnds = new int[pairs], readEnds = new int[pairs];
		SwKernel.scorePairs(params, matrix.scores(), queries.residues, queries.offsets, targets.residues, targets.offsets, h,
				e, scores, refEnds, readEnds);
		if (!Arrays.equals(scores, scorer.scores) || !Arrays.equals(refEnds, scorer.refEnds)
				|| !Arrays.equals(readEnds, scorer.readEnds)) {
			throw new IllegalStateException("device results differ from the plain Java row kernel");
		}
	}
}
