package swiss.sib.cssw4j.tornado;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Encoded sequences laid out one after another with an offset table, the form the kernel reads. */
public final class Sequences {
	public final List<String> ids;
	public final byte[] residues;
	/** Start of every sequence, plus the total length at the end. */
	public final int[] offsets;
	public final int maxLength;

	private Sequences(List<String> ids, byte[] residues, int[] offsets, int maxLength) {
		this.ids = ids;
		this.residues = residues;
		this.offsets = offsets;
		this.maxLength = maxLength;
	}

	public int count() {
		return ids.size();
	}

	/** The residues widened to ints, for kernels that must not do arithmetic on byte values (see SwKernel). */
	public int[] residuesAsInts() {
		int[] out = new int[residues.length];
		for (int i = 0; i < out.length; i++) {
			out[i] = residues[i];
		}
		return out;
	}

	/** Sequences given in memory, encoded with the matrix's letter table. */
	public static Sequences of(Matrix matrix, List<String> ids, List<String> sequences) {
		List<StringBuilder> seqs = new ArrayList<>();
		for (String s : sequences) {
			seqs.add(new StringBuilder(s));
		}
		return pack(matrix, new ArrayList<>(ids), seqs);
	}

	/** Reads FASTA files; the id of a sequence is the first word of its header. */
	public static Sequences read(Matrix matrix, List<Path> files) throws IOException {
		List<String> ids = new ArrayList<>();
		List<StringBuilder> seqs = new ArrayList<>();
		for (Path file : files) {
			try (BufferedReader in = Files.newBufferedReader(file)) {
				for (String line = in.readLine(); line != null; line = in.readLine()) {
					if (line.startsWith(">")) {
						ids.add(line.substring(1).split("[ \t]", 2)[0]);
						seqs.add(new StringBuilder());
					} else {
						seqs.get(seqs.size() - 1).append(line.strip());
					}
				}
			}
		}
		return pack(matrix, ids, seqs);
	}

	private static Sequences pack(Matrix matrix, List<String> ids, List<StringBuilder> seqs) {
		int[] offsets = new int[seqs.size() + 1];
		int total = 0, max = 0;
		for (int i = 0; i < seqs.size(); i++) {
			offsets[i] = total;
			total += seqs.get(i).length();
			max = Math.max(max, seqs.get(i).length());
		}
		offsets[seqs.size()] = total;
		byte[] residues = new byte[total];
		for (int i = 0; i < seqs.size(); i++) {
			byte[] encoded = matrix.encode(seqs.get(i).toString());
			System.arraycopy(encoded, 0, residues, offsets[i], encoded.length);
		}
		return new Sequences(ids, residues, offsets, max);
	}
}
