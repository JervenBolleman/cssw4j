package swiss.sib.cssw4j.tornado;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;

/**
 * Experiment 6: how much of the strip kernel's time goes to the substitution-matrix lookups?
 * <p>
 * Every cell reads {@code matrix[targetLetter * n + queryLetter]} from global memory. Neighbouring work-items have
 * different target letters, so a warp reads up to 32 scattered matrix entries per cell. {@link #matrixLookup} is the
 * strip kernel; {@link #noMatrixLookup} is {@link DiagnosticKernels#strip4LetterScore}, the same kernel scoring by
 * letter comparison, without the matrix reads. Its scores are different by design; the setup checks it against itself
 * run as plain Java. Balanced workload (128 x 128 random proteins of 300 residues); throughput in the secondary result
 * {@code :cells} (ops/ms; divided by 10^6 it is GCUPS).
 */
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
@State(Scope.Benchmark)
public class MatrixGatherBench {
	private TornadoScorer scorer;
	private TornadoExecutionPlan diagnostic;
	private int[] scores, refEnds, readEnds;
	private long cells;

	@Setup(Level.Trial)
	public void setup() throws Exception {
		Matrix matrix = Workloads.blosum62();
		Sequences queries = Workloads.random(matrix, "q", 128, 300, 11);
		Sequences targets = Workloads.random(matrix, "t", 128, 300, 12);
		cells = Workloads.cells(queries, targets);
		scorer = new TornadoScorer(matrix, queries, targets, Workloads.GAP_OPEN, Workloads.GAP_EXTEND,
				Layout.INTERLEAVED, TornadoScorer.Variant.STRIP4);
		scorer.run();
		Reference.check(matrix, queries, targets, scorer);

		int pairs = queries.count() * targets.count();
		int[] steps = Layout.INTERLEAVED.steps(pairs, targets.maxLength);
		int[] params = { Workloads.GAP_OPEN, Workloads.GAP_EXTEND, matrix.size(), steps[0], steps[1] };
		int[] queryInts = queries.residuesAsInts();
		int[] h = new int[(targets.maxLength + 1) * pairs], f = new int[h.length];
		scores = new int[pairs];
		refEnds = new int[pairs];
		readEnds = new int[pairs];
		TaskGraph graph = new TaskGraph("diag") //
				.transferToDevice(DataTransferMode.FIRST_EXECUTION, params, matrix.scores(), queryInts, queries.offsets,
						targets.residues, targets.offsets, h, f) //
				.task("score", DiagnosticKernels::strip4LetterScore, params, matrix.scores(), queryInts, queries.offsets,
						targets.residues, targets.offsets, h, f, scores, refEnds, readEnds) //
				.transferToHost(DataTransferMode.EVERY_EXECUTION, scores, refEnds, readEnds);
		diagnostic = new TornadoExecutionPlan(graph.snapshot()).withDevice(TornadoScorer.device());
		diagnostic.execute();

		int[] js = new int[pairs], jr = new int[pairs], jq = new int[pairs];
		DiagnosticKernels.strip4LetterScore(params, matrix.scores(), queryInts, queries.offsets, targets.residues,
				targets.offsets, new int[h.length], new int[h.length], js, jr, jq);
		if (!Arrays.equals(js, scores) || !Arrays.equals(jr, refEnds) || !Arrays.equals(jq, readEnds)) {
			throw new IllegalStateException("diagnostic kernel differs between device and plain Java");
		}
	}

	@TearDown(Level.Trial)
	public void tearDown() {
		scorer.close();
		diagnostic.freeDeviceMemory();
	}

	@Benchmark
	public void matrixLookup(CellCounter counter, Blackhole bh) {
		scorer.run();
		counter.cells += cells;
		bh.consume(scorer.scores);
	}

	@Benchmark
	public void noMatrixLookup(CellCounter counter, Blackhole bh) {
		diagnostic.execute();
		counter.cells += cells;
		bh.consume(scores);
	}
}
