package swiss.sib.cssw4j;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/**
 * A square substitution matrix together with the table that encodes sequence letters into matrix indices.
 * <p>
 * Instances are immutable and thread safe.
 */
public final class ScoringMatrix {
	/** Amino acid order of the built-in BLOSUM50 matrix (same as the SSW C library). */
	public static final String AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYVBZX*";

	/** Letter to index table for amino acids used by the SSW C library; unknown letters map to '*'. */
	private static final byte[] AA_TABLE = {
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
	};

	/** Letter to index table for nucleotides used by the SSW C library; A C G T/U = 0..3, anything else 4. */
	private static final byte[] NT_TABLE = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};

	/** BLOSUM50, as built into the SSW C library ({@code ssw_test -p}). */
	public static final ScoringMatrix BLOSUM50 = new ScoringMatrix(24, new byte[] {
		//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
		  5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5, // A
		 -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5, // R
		 -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5, // N
		 -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5, // D
		 -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5, // C
		 -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5, // Q
		 -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5, // E
		  0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5, // G
		 -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5, // H
		 -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5, // I
		 -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5, // L
		 -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5, // K
		 -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5, // M
		 -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5, // F
		 -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5, // P
		  1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5, // S
		  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, // T
		 -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, // W
		 -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, // Y
		  0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, // V
		 -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, // B
		 -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, // Z
		 -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, // X
		 -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1  // *
	}, AA_TABLE);

	private static final class Blosum62Holder {
		static final ScoringMatrix BLOSUM62 = resource("blosum62.txt");
	}

	private final int n;
	private final byte[] matrix;
	private final byte[] table;
	private final byte[] transposed;

	private ScoringMatrix(int n, byte[] matrix, byte[] table) {
		if (matrix.length < n * n) {
			throw new IllegalArgumentException("matrix needs " + n * n + " entries, got " + matrix.length);
		}
		this.n = n;
		this.matrix = Arrays.copyOf(matrix, n * n);
		this.table = table.clone();
		this.transposed = new byte[n * n];
		for (int r = 0; r < n; r++) {
			for (int c = 0; c < n; c++) {
				transposed[c * n + r] = matrix[r * n + c];
			}
		}
	}

	/** BLOSUM62 in NCBI format (the {@code blosum62.txt} distributed with the SSW C library). */
	public static ScoringMatrix blosum62() {
		return Blosum62Holder.BLOSUM62;
	}

	/**
	 * The nucleotide matrix used by {@code ssw_test}: A C G T and N (index 4) which scores 0 against everything.
	 */
	public static ScoringMatrix nucleotide(int match, int mismatch) {
		byte[] m = new byte[25];
		int k = 0;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				m[k++] = (byte) (i == j ? match : -mismatch);
			}
			m[k++] = 0;
		}
		return new ScoringMatrix(5, m, NT_TABLE);
	}

	/**
	 * Any square matrix with a caller supplied letter table.
	 *
	 * @param n      edge length of the matrix
	 * @param matrix row major {@code n * n} scores, row = reference letter, column = query letter
	 * @param table  128 entries mapping ASCII letters to indices in {@code [0, n)}
	 */
	public static ScoringMatrix of(int n, byte[] matrix, byte[] table) {
		if (table.length != 128) {
			throw new IllegalArgumentException("table must have 128 entries");
		}
		for (byte b : table) {
			if (b < 0 || b >= n) {
				throw new IllegalArgumentException("table entry " + b + " outside matrix of size " + n);
			}
		}
		return new ScoringMatrix(n, matrix, table);
	}

	/**
	 * Parses an NCBI style matrix (BLOSUM/PAM files), using the same rules as {@code ssw_test -a}: every line starting
	 * with an upper case letter or '*' is a row, and the letter is mapped to that row's index. Letters not present in
	 * the file keep their default amino acid index.
	 */
	public static ScoringMatrix parse(Reader reader) throws IOException {
		byte[] table = AA_TABLE.clone();
		byte[] m = new byte[1024];
		int k = 0, rows = 0;
		BufferedReader in = reader instanceof BufferedReader b ? b : new BufferedReader(reader);
		for (String line = in.readLine(); line != null; line = in.readLine()) {
			if (line.isEmpty()) {
				continue;
			}
			char first = line.charAt(0);
			if (first != '*' && (first < 'A' || first > 'Z')) {
				continue;
			}
			if (first != '*') {
				table[first] = table[first + 32] = (byte) rows;
			}
			int start = -1;
			for (int l = 1; l <= line.length(); ++l) {
				char c = l < line.length() ? line.charAt(l) : ' ';
				boolean digit = (c >= '0' && c <= '9') || c == '-';
				if (digit && start < 0) {
					start = l;
				} else if (!digit && start >= 0) {
					m[k++] = (byte) Integer.parseInt(line, start, l, 10);
					start = -1;
				}
			}
			++rows;
		}
		if (k == 0) {
			throw new IOException("no matrix rows found");
		}
		return new ScoringMatrix(rows, m, table);
	}

	private static ScoringMatrix resource(String name) {
		try (InputStream is = ScoringMatrix.class.getResourceAsStream(name)) {
			return parse(new InputStreamReader(is, StandardCharsets.US_ASCII));
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	/** Edge length of the square matrix. */
	public int size() {
		return n;
	}

	/** Score of reference letter index {@code ref} against query letter index {@code read}. */
	public int score(int ref, int read) {
		return matrix[ref * n + read];
	}

	/** Smallest score in the matrix. */
	public int minScore() {
		int min = 0;
		for (byte b : matrix) {
			min = Math.min(min, b);
		}
		return min;
	}

	/** Encodes ASCII letters into matrix indices. Characters outside ASCII are treated as unknown. */
	public byte[] encode(CharSequence seq) {
		byte[] out = new byte[seq.length()];
		for (int i = 0; i < out.length; ++i) {
			out[i] = table[seq.charAt(i) & 127];
		}
		return out;
	}

	/** Encodes ASCII letters into matrix indices. */
	public byte[] encode(byte[] ascii, int offset, int length) {
		byte[] out = new byte[length];
		for (int i = 0; i < length; ++i) {
			out[i] = table[ascii[offset + i] & 127];
		}
		return out;
	}

	/** The transposed matrix ({@code [read * n + ref]}), shared with the kernels; never modified. */
	byte[] transposed() {
		return transposed;
	}

	/** The row major matrix, shared with the kernels; never modified. */
	byte[] matrix() {
		return matrix;
	}
}
