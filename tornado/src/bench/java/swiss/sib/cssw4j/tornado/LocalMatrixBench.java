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
 * Experiment 7: does moving the substitution matrix into local memory recover the cost found in experiment 6?
 * <p>
 * {@code STRIP4} reads the matrix from global memory; {@code STRIP4_LOCAL} first copies it into the work-group's local
 * (shared) memory, 2.3 KB for BLOSUM62. Three workloads: {@link #balanced} (128 x 128 random proteins of 300 residues),
 * {@link #dataset} (the 3,025 pairs of the test dataset) and {@link #largestPair} (its largest pair alone). Every case is
 * checked against the plain Java row kernel in setup. Throughput in the secondary result {@code :cells} (ops/ms; divided
 * by 10^6 it is GCUPS).
 */
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
public class LocalMatrixBench {

	/** A scorer for one workload and variant, checked against the reference in setup. */
	abstract static class Case {
		TornadoScorer scorer;
		long cells;

		void build(Sequences queries, Sequences targets, Matrix matrix, String variant) {
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
	public static class Balanced extends Case {
		@Param({ "STRIP4", "STRIP4_LOCAL" })
		public String variant;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			build(Workloads.random(matrix, "q", 128, 300, 11), Workloads.random(matrix, "t", 128, 300, 12), matrix,
					variant);
		}
	}

	@State(Scope.Benchmark)
	public static class Dataset extends Case {
		@Param({ "STRIP4", "STRIP4_LOCAL" })
		public String variant;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			build(Workloads.datasetQueries(matrix), Workloads.datasetTargets(matrix), matrix, variant);
		}
	}

	@State(Scope.Benchmark)
	public static class LargestPair extends Case {
		@Param({ "STRIP4", "STRIP4_LOCAL" })
		public String variant;

		@Setup(Level.Trial)
		public void setup() throws Exception {
			Matrix matrix = Workloads.blosum62();
			Sequences[] pair = Workloads.largestPair(matrix);
			build(pair[0], pair[1], matrix, variant);
		}
	}

	@Benchmark
	public void balanced(Balanced state, CellCounter counter, Blackhole bh) {
		run(state, counter, bh);
	}

	@Benchmark
	public void dataset(Dataset state, CellCounter counter, Blackhole bh) {
		run(state, counter, bh);
	}

	@Benchmark
	public void largestPair(LargestPair state, CellCounter counter, Blackhole bh) {
		run(state, counter, bh);
	}

	private static void run(Case state, CellCounter counter, Blackhole bh) {
		state.scorer.run();
		counter.cells += state.cells;
		bh.consume(state.scorer.scores);
	}
}
