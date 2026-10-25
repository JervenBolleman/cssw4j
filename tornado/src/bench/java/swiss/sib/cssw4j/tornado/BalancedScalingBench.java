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
 * Experiment 3: how fast is the kernel when the work is balanced and there are enough pairs?
 * <p>
 * {@code n} random queries against {@code n} random targets, all 300 residues long, so every pair has 90,000 cells and
 * no pair holds up the others. {@code n} = 32, 64, 128, 256 gives 1,024 to 65,536 pairs. Throughput in the secondary
 * result {@code :cells} (ops/ms; divided by 10^6 it is GCUPS).
 */
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
@State(Scope.Benchmark)
public class BalancedScalingBench {
	static final int LENGTH = 300;

	@Param({ "32", "64", "128", "256" })
	public int n;

	private TornadoScorer scorer;
	private long cellsPerCall;

	@Setup(Level.Trial)
	public void setup() throws Exception {
		Matrix matrix = Workloads.blosum62();
		Sequences queries = Workloads.random(matrix, "q", n, LENGTH, 11);
		Sequences targets = Workloads.random(matrix, "t", n, LENGTH, 12);
		cellsPerCall = Workloads.cells(queries, targets);
		scorer = new TornadoScorer(matrix, queries, targets, Workloads.GAP_OPEN, Workloads.GAP_EXTEND, Layout.INTERLEAVED);
		scorer.run();
		Reference.check(matrix, queries, targets, scorer);
	}

	@TearDown(Level.Trial)
	public void tearDown() {
		scorer.close();
	}

	@Benchmark
	public void balanced(CellCounter counter, Blackhole bh) {
		scorer.run();
		counter.cells += cellsPerCall;
		bh.consume(scorer.scores);
	}
}
