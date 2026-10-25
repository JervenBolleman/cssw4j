package swiss.sib.cssw4j.tornado;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

/**
 * An NCBI style substitution matrix (the BLOSUM and PAM files) with the letter encoding of the SSW C library: the
 * header row gives the index of every letter, and letters that are not in the matrix map to {@code '*'}.
 */
public record Matrix(int size, int[] scores, byte[] table) {

	public static Matrix read(Path file) throws IOException {
		try (BufferedReader in = Files.newBufferedReader(file)) {
			StringBuilder letters = new StringBuilder();
			int[] scores = new int[1024];
			int k = 0;
			for (String line = in.readLine(); line != null; line = in.readLine()) {
				if (line.isEmpty()) {
					continue;
				}
				char first = line.charAt(0);
				if (first != '*' && (first < 'A' || first > 'Z')) {
					continue; // comment or header row
				}
				letters.append(first);
				for (String field : line.substring(1).trim().split("\\s+")) {
					if (!field.isEmpty()) {
						scores[k++] = Integer.parseInt(field);
					}
				}
			}
			int rows = letters.length();
			int star = letters.indexOf("*");
			byte[] table = new byte[128];
			Arrays.fill(table, (byte) Math.max(star, 0)); // unknown letters score as '*'
			for (int r = 0; r < rows; r++) {
				char c = letters.charAt(r);
				if (c != '*') {
					table[c] = table[Character.toLowerCase(c)] = (byte) r;
				}
			}
			return new Matrix(rows, Arrays.copyOf(scores, rows * rows), table);
		}
	}

	public byte[] encode(String sequence) {
		byte[] out = new byte[sequence.length()];
		for (int i = 0; i < out.length; i++) {
			out[i] = table[sequence.charAt(i) & 127];
		}
		return out;
	}
}
