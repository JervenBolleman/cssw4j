package swiss.sib.cssw4j;

import static jdk.incubator.vector.VectorOperators.AND;
import static jdk.incubator.vector.VectorOperators.LSHL;
import static jdk.incubator.vector.VectorOperators.LSHR;
import static jdk.incubator.vector.VectorOperators.OR;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorShape;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

/**
 * The vector species used by the striped kernels.
 * <p>
 * Defaults to the widest shape the CPU supports (32 byte lanes on AVX2, 64 on AVX-512). The C library is fixed at
 * 128 bits (16 byte lanes); set the system property {@code cssw4j.vectorBits=128} to use the same lane count. The
 * species are static final so the JIT can treat them as constants.
 */
final class Lanes {
	static final VectorShape SHAPE = shape();
	static final VectorSpecies<Byte> BYTE = VectorSpecies.of(byte.class, SHAPE);
	static final VectorSpecies<Short> SHORT = VectorSpecies.of(short.class, SHAPE);
	private static final VectorSpecies<Integer> INT = VectorSpecies.of(int.class, SHAPE);

	/*
	 * Shifting all lanes up by one (the _mm_slli_si128 of the C code) is done on 32 bit lanes: shift each int and
	 * carry in the top lane of the previous int (one vpermd). A byte slice()/unslice() compiles to a generic rearrange
	 * of about 30 instructions on AVX2, this to about 7.
	 */
	private static final VectorShuffle<Integer> PREVIOUS = VectorShuffle.fromOp(INT, i -> Math.max(i - 1, 0));
	private static final int[] NOT_FIRST = notFirst();
	/** Excess-128 zero (0x80) in byte lane 0, the fill of {@link #shiftBytes}. */
	private static final int[] BYTE_FILL = byteFill();

	private Lanes() {
	}

	/**
	 * Lane {@code k} of the result is lane {@code k - 1} of {@code v}; lane 0 becomes {@code Byte.MIN_VALUE}, the zero
	 * of the excess-128 byte kernel.
	 */
	static ByteVector shiftBytes(ByteVector v) {
		IntVector l = (IntVector) v.reinterpretAsInts();
		IntVector carry = l.rearrange(PREVIOUS).lanewise(LSHR, 24)
				.lanewise(AND, IntVector.fromArray(INT, NOT_FIRST, 0));
		return l.lanewise(LSHL, 8).lanewise(OR, carry).lanewise(OR, IntVector.fromArray(INT, BYTE_FILL, 0))
				.reinterpretAsBytes();
	}

	/** Lane {@code k} of the result is lane {@code k - 1} of {@code v}; lane 0 becomes 0. */
	static ShortVector shiftShorts(ShortVector v) {
		IntVector l = (IntVector) v.reinterpretAsInts();
		IntVector carry = l.rearrange(PREVIOUS).lanewise(LSHR, 16)
				.lanewise(AND, IntVector.fromArray(INT, NOT_FIRST, 0));
		return l.lanewise(LSHL, 16).lanewise(OR, carry).reinterpretAsShorts();
	}

	private static int[] notFirst() {
		int[] m = new int[INT.length()];
		java.util.Arrays.fill(m, 1, m.length, -1);
		return m;
	}

	private static int[] byteFill() {
		int[] m = new int[INT.length()];
		m[0] = 0x80;
		return m;
	}

	private static VectorShape shape() {
		String bits = System.getProperty("cssw4j.vectorBits");
		if (bits != null) {
			return VectorShape.forBitSize(Integer.parseInt(bits));
		}
		VectorShape preferred = VectorShape.preferredShape();
		// The striped algorithm needs at least two lanes of 16 bits to be worthwhile.
		return preferred.vectorBitSize() < 128 ? VectorShape.S_128_BIT : preferred;
	}
}
