package swiss.sib.cssw4j;

import java.util.Arrays;

/**
 * CIGAR helpers. Operations are packed as in BAM: the length in the high 28 bits, the operation code
 * ({@code MIDNSHP=X} = 0..8) in the low 4 bits.
 */
public final class Cigar {
	public static final String OPS = "MIDNSHP=X";
	public static final int M = 0, I = 1, D = 2, N = 3, S = 4, H = 5, P = 6, EQ = 7, X = 8;
	private static final int SHIFT = 4;

	private Cigar() {
	}

	public static int toInt(int length, int op) {
		return (length << SHIFT) | op;
	}

	public static int length(int cigarInt) {
		return cigarInt >>> SHIFT;
	}

	/** The operation code, {@link #M} for codes outside {@code 0..8} (as in the C library). */
	public static int op(int cigarInt) {
		int op = cigarInt & 0xf;
		return op > 8 ? M : op;
	}

	public static char opChar(int cigarInt) {
		return OPS.charAt(op(cigarInt));
	}

	public static String toString(int[] cigar) {
		if (cigar == null || cigar.length == 0) {
			return "*";
		}
		StringBuilder sb = new StringBuilder(cigar.length * 4);
		for (int c : cigar) {
			sb.append(length(c)).append(opChar(c));
		}
		return sb.toString();
	}

	/** Result of {@link #markMismatch}: the extended cigar and the edit distance. */
	public record Marked(int mismatches, int[] cigar) {
		@Override
		public String toString() {
			return Cigar.toString(cigar) + " NM:" + mismatches;
		}
	}

	/**
	 * Port of {@code mark_mismatch}: rewrites {@code M} into {@code =}/{@code X}, adds soft clips and counts the
	 * mismatches plus inserted and deleted letters.
	 *
	 * @param ref     encoded reference
	 * @param read    encoded query
	 * @param readLen length of the full query
	 */
	public static Marked markMismatch(Alignment a, byte[] ref, byte[] read, int readLen) {
		return markMismatch(a.refBegin1(), a.readBegin1(), a.readEnd1(), ref, read, readLen, a.cigar());
	}

	public static Marked markMismatch(int refBegin1, int readBegin1, int readEnd1, byte[] ref, byte[] read,
			int readLen, int[] cigar) {
		int[] out = new int[cigar.length + 2];
		int p = 0, mismatches = 0, lengthM = 0, lengthX = 0;
		int r = refBegin1, q = readBegin1;
		if (readBegin1 > 0) {
			out[p++] = toInt(readBegin1, S);
		}
		for (int c : cigar) {
			int op = op(c), len = length(c);
			if (op == M) {
				for (int j = 0; j < len; ++j, ++r, ++q) {
					if (ref[r] != read[q]) {
						++mismatches;
						if (lengthM != 0) {
							out = add(out, p++, lengthM, EQ);
							lengthM = 0;
						}
						++lengthX;
					} else {
						if (lengthX != 0) {
							out = add(out, p++, lengthX, X);
							lengthX = 0;
						}
						++lengthM;
					}
				}
			} else if (op == I || op == D) {
				if (op == I) {
					q += len;
				} else {
					r += len;
				}
				mismatches += len;
				if (lengthM != 0) {
					out = add(out, p++, lengthM, EQ);
					lengthM = 0;
				} else if (lengthX != 0) {
					out = add(out, p++, lengthX, X);
					lengthX = 0;
				}
				out = add(out, p++, len, op);
			}
		}
		if (lengthM != 0) {
			out = add(out, p++, lengthM, EQ);
		} else if (lengthX != 0) {
			out = add(out, p++, lengthX, X);
		}
		int clip = readLen - readEnd1 - 1;
		if (clip > 0) {
			out = add(out, p++, clip, S);
		}
		return new Marked(mismatches, Arrays.copyOf(out, p));
	}

	private static int[] add(int[] out, int p, int length, int op) {
		if (p >= out.length) {
			out = Arrays.copyOf(out, out.length * 2);
		}
		out[p] = toInt(length, op);
		return out;
	}
}
