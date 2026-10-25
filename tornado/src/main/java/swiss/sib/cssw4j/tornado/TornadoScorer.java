package swiss.sib.cssw4j.tornado;

import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;

/**
 * Runs {@link SwKernel#scorePairs} on a TornadoVM device for every pair of a batch of queries and a batch of targets.
 * The task graph is built once; {@link #run()} can be called repeatedly (the first call compiles the kernel).
 */
public final class TornadoScorer implements AutoCloseable {
	public final int[] scores, refEnds, readEnds;
	private final TornadoExecutionPlan plan;

	public TornadoScorer(Matrix matrix, Sequences queries, Sequences targets, int gapOpen, int gapExtend,
			Layout layout) {
		int pairs = queries.count() * targets.count();
		int[] steps = layout.steps(pairs, queries.maxLength);
		int[] params = { gapOpen, gapExtend, matrix.size(), steps[0], steps[1] };
		// Scratch for the recurrence: one H and one E value per query position and pair, interleaved by pair.
		int[] h = new int[(queries.maxLength + 1) * pairs];
		int[] e = new int[(queries.maxLength + 1) * pairs];
		scores = new int[pairs];
		refEnds = new int[pairs];
		readEnds = new int[pairs];

		TaskGraph graph = new TaskGraph("sw") //
				.transferToDevice(DataTransferMode.FIRST_EXECUTION, params, matrix.scores(), queries.residues,
						queries.offsets, targets.residues, targets.offsets, h, e) //
				.task("score", SwKernel::scorePairs, params, matrix.scores(), queries.residues, queries.offsets,
						targets.residues, targets.offsets, h, e, scores, refEnds, readEnds) //
				.transferToHost(DataTransferMode.EVERY_EXECUTION, scores, refEnds, readEnds);
		ImmutableTaskGraph immutable = graph.snapshot();
		plan = new TornadoExecutionPlan(immutable).withDevice(device());
	}

	/**
	 * The device named by the system property {@code cssw4j.device} as {@code driver:index} (for example {@code 0:1}),
	 * from {@code tornado --devices}. Without the property the TornadoVM default device is used.
	 */
	public static uk.ac.manchester.tornado.api.common.TornadoDevice device() {
		String spec = System.getProperty("cssw4j.device", "0:0");
		String[] parts = spec.split(":");
		return TornadoExecutionPlan.getDevice(Integer.parseInt(parts[0]), Integer.parseInt(parts[1]));
	}

	public void run() {
		plan.execute();
	}

	@Override
	public void close() {
		plan.freeDeviceMemory();
	}
}
