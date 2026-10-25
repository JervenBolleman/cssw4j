package swiss.sib.cssw4j.tornado;

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
 * Experiment 1: is the GPU run slow because of host-device transfers or kernel launch overhead?
 * <p>
 * {@link #fullRun} is the scoring task graph on the test dataset (3,025 pairs). {@link #overheadOnly} is a task graph
 * with the same arguments and transfers whose kernel only writes the result arrays, so it measures dispatch, transfers
 * and the copy back without the alignment work. Average time per call.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
@State(Scope.Benchmark)
public class TransferBench {
	private TornadoScorer scorer;
	private TornadoExecutionPlan overhead;
	private int[] scores, refEnds, readEnds;

	@Setup(Level.Trial)
	public void setup() throws Exception {
		Matrix matrix = Workloads.blosum62();
		Sequences queries = Workloads.datasetQueries(matrix);
		Sequences targets = Workloads.datasetTargets(matrix);
		scorer = new TornadoScorer(matrix, queries, targets, Workloads.GAP_OPEN, Workloads.GAP_EXTEND, Layout.INTERLEAVED);
		scorer.run();
		Reference.check(matrix, queries, targets, scorer);

		int pairs = queries.count() * targets.count();
		int[] steps = Layout.INTERLEAVED.steps(pairs, queries.maxLength);
		int[] params = { Workloads.GAP_OPEN, Workloads.GAP_EXTEND, matrix.size(), steps[0], steps[1] };
		int[] h = new int[(queries.maxLength + 1) * pairs], e = new int[h.length];
		scores = new int[pairs];
		refEnds = new int[pairs];
		readEnds = new int[pairs];
		TaskGraph graph = new TaskGraph("overhead") //
				.transferToDevice(DataTransferMode.FIRST_EXECUTION, params, matrix.scores(), queries.residues,
						queries.offsets, targets.residues, targets.offsets, h, e) //
				.task("clear", OverheadKernel::clear, params, matrix.scores(), queries.residues, queries.offsets,
						targets.residues, targets.offsets, h, e, scores, refEnds, readEnds) //
				.transferToHost(DataTransferMode.EVERY_EXECUTION, scores, refEnds, readEnds);
		overhead = new TornadoExecutionPlan(graph.snapshot()).withDevice(TornadoScorer.device());
		overhead.execute();
	}

	@TearDown(Level.Trial)
	public void tearDown() {
		scorer.close();
		overhead.freeDeviceMemory();
	}

	@Benchmark
	public void fullRun(Blackhole bh) {
		scorer.run();
		bh.consume(scorer.scores);
	}

	@Benchmark
	public void overheadOnly(Blackhole bh) {
		overhead.execute();
		bh.consume(scores);
	}
}
