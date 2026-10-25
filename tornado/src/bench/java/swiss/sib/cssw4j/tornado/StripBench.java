package swiss.sib.cssw4j.tornado;

import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Experiment 5: is the kernel limited by the memory traffic of its H and E rows?
 * <p>
 * The row kernel reads and writes two ints per cell in global memory. The strip kernel keeps four query rows in
 * registers and writes only the bottom row of each strip, a quarter of that traffic. If the row kernel were bound by
 * memory bandwidth, the strip kernel would be up to four times faster.
 * <p>
 * {@link #balanced}: 128 x 128 random proteins of 300 residues (experiment 3). {@link #largestPair}: the largest pair
 * of the test dataset alone (experiment 2), which shows the speed of one GPU thread. Throughput in the secondary result
 * {@code :cells} (ops/ms; divided by 10^6 it is GCUPS).
 */
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
public class StripBench {

	@State(Scope.Benchmark)
	public static class Balanced {
		@Param({ "ROW", "STRIP4" })
		public String variant;
		TornadoScorer scorer;
		long cells;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			Sequences queries = Workloads.random(matrix, "q", 128, 300, 11);
			Sequences targets = Workloads.random(matrix, "t", 128, 300, 12);
			cells = Workloads.cells(queries, targets);
			scorer = new TornadoScorer(matrix, queries, targets, Workloads.GAP_OPEN, Workloads.GAP_EXTEND,
					Layout.INTERLEAVED, TornadoScorer.Variant.parse(variant));
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
		@Param({ "ROW", "STRIP4" })
		public String variant;
		TornadoScorer scorer;
		long cells;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			Sequences[] pair = Workloads.largestPair(matrix);
			cells = Workloads.cells(pair[0], pair[1]);
			scorer = new TornadoScorer(matrix, pair[0], pair[1], Workloads.GAP_OPEN, Workloads.GAP_EXTEND,
					Layout.INTERLEAVED, TornadoScorer.Variant.parse(variant));
			scorer.run();
			Reference.check(matrix, pair[0], pair[1], scorer);
		}

		@TearDown(Level.Trial)
		public void tearDown() {
			scorer.close();
		}
	}

	@Benchmark
	public void balanced(Balanced state, CellCounter counter, Blackhole bh) {
		state.scorer.run();
		counter.cells += state.cells;
		bh.consume(state.scorer.scores);
	}

	@Benchmark
	public void largestPair(LargestPair state, CellCounter counter, Blackhole bh) {
		state.scorer.run();
		counter.cells += state.cells;
		bh.consume(state.scorer.scores);
	}
}
