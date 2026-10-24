package swiss.sib.cssw4j;

import java.util.Arrays;
import java.util.Random;
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
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Strategy: encode the traceback direction of each cell without branches.
 * <p>
 * Why: for every cell of the banded alignment the traceback needs to know which neighbour the best score came from.
 * The direct way is a chain of {@code if}/{@code else} on the comparison results. For real sequences those
 * comparisons are close to random, so the CPU mispredicts about half of the branches. The branch-free version turns
 * each comparison into a sign bit ({@code (a - b) >>> 31}) and builds the direction with arithmetic.
 * <p>
 * Result in the port: full-CIGAR throughput went from about 1.5 to 1.9 GCUPS after this change (the number includes
 * the other row-loop changes made at the same time, see the commit history).
 * <p>
 * The inputs are random and rotate through {@link #BLOCKS} different blocks. Repeating one block would let the branch
 * predictor learn it, which hides the cost of mispredictions. The setup checks that both versions produce the same
 * bytes.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(2)
@State(Scope.Thread)
public class DirectionEncodingBench {
	/** Same bit layout as {@link Aligner}. */
	static final int DIR_E_OPEN = 8, DIR_F_OPEN = 16;

	/** Number of distinct input blocks; each call uses the next one, so the branch predictor cannot memorize them. */
	static final int BLOCKS = 16;

	@Param({ "256", "4096" })
	int cells;

	int[] eA, eB, fA, fB, e1, f1, diag;
	byte[] outBranchy, outBranchless;
	int call;

	@Setup(Level.Trial)
	public void setup() {
		Random rnd = new Random(42);
		eA = random(rnd, cells * BLOCKS);
		eB = random(rnd, cells * BLOCKS);
		fA = random(rnd, cells * BLOCKS);
		fB = random(rnd, cells * BLOCKS);
		e1 = random(rnd, cells * BLOCKS);
		f1 = random(rnd, cells * BLOCKS);
		diag = random(rnd, cells * BLOCKS);
		outBranchy = new byte[cells];
		outBranchless = new byte[cells];
		branchy(eA, eB, fA, fB, e1, f1, diag, 0, outBranchy);
		branchless(eA, eB, fA, fB, e1, f1, diag, 0, outBranchless);
		if (!Arrays.equals(outBranchy, outBranchless)) {
			throw new IllegalStateException("the two encodings differ");
		}
	}

	/** Start of the next input block. */
	int nextOffset() {
		return (call++ % BLOCKS) * cells;
	}

	static int[] random(Random rnd, int n) {
		int[] a = new int[n];
		for (int i = 0; i < n; i++) {
			a[i] = rnd.nextInt(101) - 50;
		}
		return a;
	}

	/** The direct version: one comparison per branch. */
	static void branchy(int[] eA, int[] eB, int[] fA, int[] fB, int[] e1, int[] f1, int[] diag, int off,
			byte[] out) {
		for (int i = 0; i < out.length; i++) {
			int k = off + i;
			int dirE = eA[k] > eB[k] ? DIR_E_OPEN : 0;
			int dirF = fA[k] > fB[k] ? DIR_F_OPEN : 0;
			int temp1 = Math.max(e1[k], f1[k]);
			int temp2 = diag[k];
			int dirH;
			if (temp1 <= temp2) {
				dirH = 1;
			} else if (e1[k] > f1[k]) {
				dirH = dirE != 0 ? 3 : 2;
			} else {
				dirH = dirF != 0 ? 5 : 4;
			}
			out[i] = (byte) (dirH | dirE | dirF);
		}
	}

	/** The version used in {@link Aligner#bandRow}: sign bits instead of branches. */
	static void branchless(int[] eA, int[] eB, int[] fA, int[] fB, int[] e1, int[] f1, int[] diag, int off,
			byte[] out) {
		for (int i = 0; i < out.length; i++) {
			int k = off + i;
			final int eOpen = (eB[k] - eA[k]) >>> 31; // 1 when eA > eB
			final int fOpen = (fB[k] - fA[k]) >>> 31;
			final int fromE = (f1[k] - e1[k]) >>> 31; // 1 when e1 > f1
			final int temp1 = Math.max(e1[k], f1[k]);
			final int diagonal = ((diag[k] - temp1) >>> 31) ^ 1; // 1 when temp1 <= diag
			final int gap = 4 - 2 * fromE + fOpen + fromE * (eOpen - fOpen);
			final int dirH = gap + diagonal * (1 - gap);
			out[i] = (byte) (dirH | eOpen << 3 | fOpen << 4);
		}
	}

	@Benchmark
	public void directIfElse(Blackhole bh) {
		branchy(eA, eB, fA, fB, e1, f1, diag, nextOffset(), outBranchy);
		bh.consume(outBranchy);
	}

	@Benchmark
	public void branchFree(Blackhole bh) {
		branchless(eA, eB, fA, fB, e1, f1, diag, nextOffset(), outBranchless);
		bh.consume(outBranchless);
	}
}
