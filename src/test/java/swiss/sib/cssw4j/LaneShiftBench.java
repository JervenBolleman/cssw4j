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
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

/**
 * Strategy: shift a byte vector up by one lane, the Vector API replacement for {@code _mm_slli_si128}, in four ways.
 * <p>
 * Why: the striped kernel shifts the vector of the last segment into lane 0 once per column, and the lazy-F loop
 * shifts {@code vF} once per round, so the shift runs very often. There is no direct shift in the Vector API, so the
 * variants tried were:
 * <ol>
 * <li>{@code slice}: the obvious call, {@code zero.slice(n - 1, v)}. Generic rearrangement, about 30 instructions on
 * AVX2.</li>
 * <li>{@code unslice + blend}: {@code v.unslice(1)} with the first lane blended from the fill. Also generic.</li>
 * <li>64 bit lanes: shift each long and carry the top byte of the previous long in with a permute.</li>
 * <li>32 bit lanes: the same on ints, with a single {@code vpermd}. This is what {@code Lanes.shiftBytes} uses.</li>
 * </ol>
 * Result: in this isolated benchmark the four variants are within about 12% of each other, and the generic slice is
 * the fastest. So the shift itself is not where the port spends its time. The port's end-to-end numbers (length-16
 * queries 0.43 to 0.60 GCUPS, length-64 queries 1.28 to 1.58 GCUPS) were taken while other changes were made too, so
 * they do not isolate this change. The 32 bit version is kept because it has the fewest instructions, not because it
 * was measured faster here.
 * <p>
 * The setup checks that every variant produces the same result, with fill value {@code Byte.MIN_VALUE} in lane 0 as
 * in the kernel.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(value = 2, jvmArgs = { "--add-modules", "jdk.incubator.vector" })
@State(Scope.Thread)
public class LaneShiftBench {
	static final VectorSpecies<Byte> BYTE = VectorSpecies.of(byte.class, Lanes.SHAPE);
	static final VectorSpecies<Integer> INT = VectorSpecies.of(int.class, Lanes.SHAPE);
	static final VectorSpecies<Long> LONG = VectorSpecies.of(long.class, Lanes.SHAPE);

	static final int VECTORS = 64;
	static final ByteVector FILL = ByteVector.broadcast(BYTE, Byte.MIN_VALUE);
	static final VectorMask<Byte> FIRST_LANE = VectorMask.fromLong(BYTE, 1L);

	/** One value per lane; the same data is shifted by every variant. */
	byte[][] data;
	byte[] result;

	@Setup(Level.Trial)
	public void setup() {
		Random rnd = new Random(7);
		data = new byte[VECTORS][BYTE.length()];
		for (byte[] v : data) {
			rnd.nextBytes(v);
		}
		ByteVector a = sliceShift(ByteVector.fromArray(BYTE, data[0], 0));
		check(a, unsliceBlend(ByteVector.fromArray(BYTE, data[0], 0)), "unslice+blend");
		check(a, longShift(ByteVector.fromArray(BYTE, data[0], 0)), "64 bit lanes");
		check(a, intShift(ByteVector.fromArray(BYTE, data[0], 0)), "32 bit lanes");
		result = new byte[BYTE.length()];
	}

	static void check(ByteVector expected, ByteVector actual, String what) {
		if (!expected.eq(actual).allTrue()) {
			throw new IllegalStateException(what + " differs from slice");
		}
	}

	/** Strategy 1: the generic slice. */
	static ByteVector sliceShift(ByteVector v) {
		return FILL.slice(BYTE.length() - 1, v);
	}

	/** Strategy 2: unslice and blend the first lane. */
	static ByteVector unsliceBlend(ByteVector v) {
		return v.unslice(1).blend(FILL, FIRST_LANE);
	}

	/** Strategy 3: 64 bit lanes. */
	static final VectorShuffle<Long> PREVIOUS_LONG = VectorShuffle.fromOp(LONG, i -> Math.max(i - 1, 0));
	static final long[] NOT_FIRST_LONG = notFirstLong();
	static final long[] FILL_LONG = fillLong();

	static long[] notFirstLong() {
		long[] m = new long[LONG.length()];
		Arrays.fill(m, 1, m.length, -1L);
		return m;
	}

	static long[] fillLong() {
		long[] m = new long[LONG.length()];
		m[0] = Byte.MIN_VALUE & 0xff;
		return m;
	}

	static ByteVector longShift(ByteVector v) {
		LongVector l = (LongVector) v.reinterpretAsLongs();
		LongVector carry = l.rearrange(PREVIOUS_LONG).lanewise(VectorOperators.LSHR, 56)
				.lanewise(VectorOperators.AND, LongVector.fromArray(LONG, NOT_FIRST_LONG, 0));
		return l.lanewise(VectorOperators.LSHL, 8).lanewise(VectorOperators.OR, carry)
				.lanewise(VectorOperators.OR, LongVector.fromArray(LONG, FILL_LONG, 0)).reinterpretAsBytes();
	}

	/** Strategy 4: 32 bit lanes, as in {@code Lanes.shiftBytes}. */
	static ByteVector intShift(ByteVector v) {
		return Lanes.shiftBytes(v);
	}

	// Each benchmark shifts the same vectors in a chain, so the dependency between shifts is measured.

	@Benchmark
	public void sliceGeneric(Blackhole bh) {
		chain(bh, 0);
	}

	@Benchmark
	public void unsliceAndBlend(Blackhole bh) {
		chain(bh, 1);
	}

	@Benchmark
	public void longLanes(Blackhole bh) {
		chain(bh, 2);
	}

	@Benchmark
	public void intLanes(Blackhole bh) {
		chain(bh, 3);
	}

	private void chain(Blackhole bh, int variant) {
		ByteVector acc = ByteVector.zero(BYTE);
		for (byte[] row : data) {
			ByteVector x = ByteVector.fromArray(BYTE, row, 0);
			ByteVector shifted = switch (variant) {
				case 0 -> sliceShift(acc);
				case 1 -> unsliceBlend(acc);
				case 2 -> longShift(acc);
				default -> intShift(acc);
			};
			acc = shifted.lanewise(VectorOperators.MAX, x);
		}
		acc.intoArray(result, 0);
		bh.consume(result);
	}
}
