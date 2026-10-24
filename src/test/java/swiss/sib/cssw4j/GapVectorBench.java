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
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Strategy: get the gap penalty vectors into the striped inner loop without a broadcast inside the loop.
 * <p>
 * Why: the inner loop subtracts the gap open and gap extend penalties from every vector, and each of those values is a
 * broadcast of one number. When the kernel was built with {@code ByteVector.broadcast} inside the loop body, the
 * assembly showed a {@code vpbroadcastb} on every iteration (C2 did not hoist the broadcast out of the loop). Variants:
 * <ol>
 * <li>{@code broadcastEachUse}: call {@code broadcast} inside the loop, the first version of the kernel.</li>
 * <li>{@code broadcastHoisted}: create the vectors once, before the loop. This is the obvious fix.</li>
 * <li>{@code loadFromArray}: write the two penalties into a small array once per call and load them as vectors. This
 * is what the port uses.</li>
 * </ol>
 * Result: in this loop the three variants are within about 5% of each other, and the array load is the slowest of
 * the three. The loop is too simple to reproduce the problem: the penalties are constants, so C2 hoists the broadcast
 * out of the loop anyway. In the real kernel the broadcasts were re-created inside the loop, which we believe is C2
 * giving up its registers under pressure and rematerializing constants (not confirmed with a profile). The port keeps
 * the array load because the change removed the broadcasts from the kernel's assembly, and the long-query pair went
 * from 7.5 to 8.9 GCUPS after it. That end-to-end number is the only evidence for the change, and it is not isolated.
 * <p>
 * The loop here follows the shape of the striped inner loop: subtract a penalty from the vector, take a maximum,
 * store. The setup checks that all three variants produce the same array.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(value = 2, jvmArgs = { "--add-modules", "jdk.incubator.vector" })
@State(Scope.Thread)
public class GapVectorBench {
	static final VectorSpecies<Byte> BYTE = VectorSpecies.of(byte.class, Lanes.SHAPE);

	/** Number of vectors per column, as for a query of several thousand residues. */
	static final int SEGMENTS = 256;
	static final int GAP_OPEN = 10, GAP_EXTEND = 2;

	byte[] column, out, gapScratch;

	@Setup(Level.Trial)
	public void setup() {
		Random rnd = new Random(3);
		column = new byte[SEGMENTS * BYTE.length()];
		for (int i = 0; i < column.length; i++) {
			column[i] = (byte) (rnd.nextInt(256) - 128);
		}
		out = new byte[column.length];
		gapScratch = new byte[2 * BYTE.length()];
		byte[] a = new byte[column.length], b = new byte[column.length], c = new byte[column.length];
		eachUse(column, a);
		hoisted(column, b);
		fromArray(column, c, gapScratch);
		if (!Arrays.equals(a, b) || !Arrays.equals(a, c)) {
			throw new IllegalStateException("variants differ");
		}
	}

	static void eachUse(byte[] in, byte[] out) {
		for (int off = 0; off < in.length; off += BYTE.length()) {
			ByteVector x = ByteVector.fromArray(BYTE, in, off);
			x = x.lanewise(VectorOperators.SSUB, ByteVector.broadcast(BYTE, (byte) GAP_OPEN))
					.lanewise(VectorOperators.MAX, ByteVector.broadcast(BYTE, (byte) GAP_EXTEND));
			x.intoArray(out, off);
		}
	}

	static void hoisted(byte[] in, byte[] out) {
		ByteVector gapOpen = ByteVector.broadcast(BYTE, (byte) GAP_OPEN);
		ByteVector gapExtend = ByteVector.broadcast(BYTE, (byte) GAP_EXTEND);
		for (int off = 0; off < in.length; off += BYTE.length()) {
			ByteVector x = ByteVector.fromArray(BYTE, in, off);
			x = x.lanewise(VectorOperators.SSUB, gapOpen).lanewise(VectorOperators.MAX, gapExtend);
			x.intoArray(out, off);
		}
	}

	static void fromArray(byte[] in, byte[] out, byte[] gaps) {
		int n = BYTE.length();
		Arrays.fill(gaps, 0, n, (byte) GAP_OPEN);
		Arrays.fill(gaps, n, 2 * n, (byte) GAP_EXTEND);
		ByteVector gapOpen = ByteVector.fromArray(BYTE, gaps, 0);
		ByteVector gapExtend = ByteVector.fromArray(BYTE, gaps, n);
		for (int off = 0; off < in.length; off += n) {
			ByteVector x = ByteVector.fromArray(BYTE, in, off);
			x = x.lanewise(VectorOperators.SSUB, gapOpen).lanewise(VectorOperators.MAX, gapExtend);
			x.intoArray(out, off);
		}
	}

	@Benchmark
	public void broadcastEachUse(Blackhole bh) {
		eachUse(column, out);
		bh.consume(out);
	}

	@Benchmark
	public void broadcastHoisted(Blackhole bh) {
		hoisted(column, out);
		bh.consume(out);
	}

	@Benchmark
	public void loadFromArray(Blackhole bh) {
		fromArray(column, out, gapScratch);
		bh.consume(out);
	}
}
