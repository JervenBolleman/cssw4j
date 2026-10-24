package swiss.sib.cssw4j;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/** Minimal FASTA reader for the test dataset; ids are the first word of the header, like the C harness. */
record Fasta(String id, String seq) {

	static List<Fasta> resource(String... names) {
		List<Fasta> out = new ArrayList<>();
		for (String name : names) {
			try (InputStream is = Fasta.class.getResourceAsStream("/dataset/" + name)) {
				if (is == null) {
					throw new IllegalArgumentException("missing test resource " + name);
				}
				read(new BufferedReader(new InputStreamReader(is, StandardCharsets.US_ASCII)), out);
			} catch (IOException e) {
				throw new UncheckedIOException(e);
			}
		}
		return out;
	}

	private static void read(BufferedReader in, List<Fasta> out) throws IOException {
		String id = null;
		StringBuilder seq = new StringBuilder();
		for (String line = in.readLine(); line != null; line = in.readLine()) {
			if (line.startsWith(">")) {
				if (id != null) {
					out.add(new Fasta(id, seq.toString()));
				}
				id = line.substring(1).split("[ \t]", 2)[0];
				seq.setLength(0);
			} else {
				seq.append(line.strip());
			}
		}
		if (id != null) {
			out.add(new Fasta(id, seq.toString()));
		}
	}
}
