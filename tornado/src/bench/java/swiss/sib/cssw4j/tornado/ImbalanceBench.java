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

/**
 * Experiment 2: does one large pair set the run time of the whole dataset?
 * <p>
 * The kernel runs one work-item per pair, and every pair is computed serially by its work-item. In the test dataset the
 * largest pair has 35 times the cells of the mean pair. {@link #allPairs} runs all 3,025 pairs; {@link #largestPairOnly}
 * runs that one pair alone. If the two times are close, the run time is the time of one GPU thread working through the
 * largest pair. Average time per call.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
public class ImbalanceBench {

	@State(Scope.Benchmark)
	public static class AllPairs {
		TornadoScorer scorer;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			Sequences queries = Workloads.datasetQueries(matrix), targets = Workloads.datasetTargets(matrix);
			scorer = new TornadoScorer(matrix, queries, targets, Workloads.GAP_OPEN, Workloads.GAP_EXTEND,
					Layout.INTERLEAVED);
			scorer.run();
			Reference.check(matrix, queries, targets, scorer);
		}

		@TearDown(Level.Trial)
		public void tearDown() {
			scorer.close();
		}
	}

	@State(Scope.Benchmark)
	public static class LargestPair {
		TornadoScorer scorer;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			Sequences[] pair = Workloads.largestPair(matrix);
			scorer = new TornadoScorer(matrix, pair[0], pair[1], Workloads.GAP_OPEN, Workloads.GAP_EXTEND,
					Layout.INTERLEAVED);
			scorer.run();
			Reference.check(matrix, pair[0], pair[1], scorer);
		}

		@TearDown(Level.Trial)
		public void tearDown() {
			scorer.close();
		}
	}

	@Benchmark
	public void allPairs(AllPairs state, Blackhole bh) {
		state.scorer.run();
		bh.consume(state.scorer.scores);
	}

	@Benchmark
	public void largestPairOnly(LargestPair state, Blackhole bh) {
		state.scorer.run();
		bh.consume(state.scorer.scores);
	}
}
