package swiss.sib.cssw4j.tornado;

/**
 * Where the kernel keeps the per-pair recurrence rows (H and E, one value per query position) in the scratch arrays.
 * <ul>
 * <li>{@link #INTERLEAVED}: the values of one query position for all pairs are adjacent, so neighbouring work-items
 * read neighbouring words. This is the layout a GPU wants (coalesced loads).</li>
 * <li>{@link #CONTIGUOUS}: the values of one pair are adjacent, which is what a CPU wants (each thread walks its own
 * cache lines). It is the better choice for a CPU OpenCL device.</li>
 * </ul>
 * The results are identical; only the speed differs.
 */
public enum Layout {
	INTERLEAVED, CONTIGUOUS;

	/** The {@code jStep} and {@code pStep} kernel parameters for this layout. */
	int[] steps(int pairs, int maxQueryLength) {
		return this == INTERLEAVED ? new int[] { pairs, 1 } : new int[] { 1, maxQueryLength + 1 };
	}

	public static Layout parse(String name) {
		return valueOf(name.toUpperCase());
	}
}
