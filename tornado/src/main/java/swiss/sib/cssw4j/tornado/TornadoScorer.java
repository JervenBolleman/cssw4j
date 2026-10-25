package swiss.sib.cssw4j.tornado;

import uk.ac.manchester.tornado.api.GridScheduler;
import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.WorkerGrid;
import uk.ac.manchester.tornado.api.WorkerGrid1D;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;

/**
 * Runs the scoring kernel on a TornadoVM device for every pair of a batch of queries and a batch of targets. The task
 * graph is built once; {@link #run()} can be called repeatedly (the first call compiles the kernel).
 */
public final class TornadoScorer implements AutoCloseable {
	/** Which kernel to run; all compute the same scores and end positions. */
	public enum Variant {
		/** {@link SwKernel#scorePairs}: one row per step, H and E rows in global memory. */
		ROW,
		/** {@link SwKernel#scorePairsStrip4}: four rows per pass in registers, a quarter of the memory traffic. */
		STRIP4,
		/** {@link SwKernel#scorePairsStrip4Local}: as STRIP4, with the substitution matrix in local memory. */
		STRIP4_LOCAL;

		public static Variant parse(String name) {
			return valueOf(name.toUpperCase());
		}
	}

	/** Work-group size for {@link Variant#STRIP4_LOCAL}. */
	static final int LOCAL_SIZE = 64;

	public final int[] scores, refEnds, readEnds;
	private final TornadoExecutionPlan plan;

	public TornadoScorer(Matrix matrix, Sequences queries, Sequences targets, int gapOpen, int gapExtend,
			Layout layout) {
		this(matrix, queries, targets, gapOpen, gapExtend, layout, Variant.ROW);
	}

	public TornadoScorer(Matrix matrix, Sequences queries, Sequences targets, int gapOpen, int gapExtend,
			Layout layout, Variant variant) {
		int pairs = queries.count() * targets.count();
		// Scratch: ROW keeps one H and one E value per query position and pair, STRIP4 one H and one F value per
		// target position and pair.
		int length = variant == Variant.ROW ? queries.maxLength : targets.maxLength;
		int[] steps = layout.steps(pairs, length);
		int[] params = { gapOpen, gapExtend, matrix.size(), steps[0], steps[1] };
		int[] h = new int[(length + 1) * pairs];
		int[] e = new int[(length + 1) * pairs];
		scores = new int[pairs];
		refEnds = new int[pairs];
		readEnds = new int[pairs];

		TaskGraph graph = new TaskGraph("sw");
		if (variant == Variant.ROW) {
			graph.transferToDevice(DataTransferMode.FIRST_EXECUTION, params, matrix.scores(), queries.residues,
					queries.offsets, targets.residues, targets.offsets, h, e);
			graph.task("score", SwKernel::scorePairs, params, matrix.scores(), queries.residues, queries.offsets,
					targets.residues, targets.offsets, h, e, scores, refEnds, readEnds);
		} else {
			int[] queryInts = queries.residuesAsInts();
			graph.transferToDevice(DataTransferMode.FIRST_EXECUTION, params, matrix.scores(), queryInts,
					queries.offsets, targets.residues, targets.offsets, h, e);
			if (variant == Variant.STRIP4) {
				graph.task("score", SwKernel::scorePairsStrip4, params, matrix.scores(), queryInts, queries.offsets,
						targets.residues, targets.offsets, h, e, scores, refEnds, readEnds);
			} else {
				if (matrix.size() * matrix.size() > SwKernel.MAX_MATRIX_CELLS) {
					throw new IllegalArgumentException("matrix too large for local memory: " + matrix.size());
				}
				graph.task("score", SwKernel::scorePairsStrip4Local, new KernelContext(), params, matrix.scores(),
						queryInts, queries.offsets, targets.residues, targets.offsets, h, e, scores, refEnds, readEnds);
			}
		}
		graph.transferToHost(DataTransferMode.EVERY_EXECUTION, scores, refEnds, readEnds);
		ImmutableTaskGraph immutable = graph.snapshot();
		TornadoExecutionPlan p = new TornadoExecutionPlan(immutable).withDevice(device());
		if (variant == Variant.STRIP4_LOCAL) {
			// One work-item per pair, in work-groups of LOCAL_SIZE that share the matrix copy; the global size must be
			// a multiple of the work-group size, the kernel ignores the extra work-items.
			WorkerGrid grid = new WorkerGrid1D((pairs + LOCAL_SIZE - 1) / LOCAL_SIZE * LOCAL_SIZE);
			grid.setLocalWork(LOCAL_SIZE, 1, 1);
			p = p.withGridScheduler(new GridScheduler("sw.score", grid));
		}
		plan = p;
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
