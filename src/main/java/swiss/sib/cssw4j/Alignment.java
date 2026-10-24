package swiss.sib.cssw4j;

/**
 * Result of {@link Aligner#align}, the Java counterpart of {@code s_align}. All positions are 0-based and inclusive.
 *
 * @param score1     the best alignment score
 * @param score2     sub-optimal alignment score (0 when maskLen &lt; 15)
 * @param refBegin1  begin of the best alignment on the reference; -1 when not computed
 * @param refEnd1    end of the best alignment on the reference
 * @param readBegin1 begin of the best alignment on the query; -1 when not computed
 * @param readEnd1   end of the best alignment on the query
 * @param refEnd2    end of the sub-optimal alignment on the reference (-1 when maskLen &lt; 15)
 * @param cigar      packed cigar ({@link Cigar}) of the best alignment; {@code null} when not computed
 * @param flag       0: path accurate; 1: banded alignment failed; 2: path may miss a small part
 */
public record Alignment(int score1, int score2, int refBegin1, int refEnd1, int readBegin1, int readEnd1, int refEnd2,
		int[] cigar, int flag) {

	public String cigarString() {
		return Cigar.toString(cigar);
	}

	@Override
	public String toString() {
		return "Alignment[score1=" + score1 + ", score2=" + score2 + ", ref=" + refBegin1 + ".." + refEnd1 + ", read="
				+ readBegin1 + ".." + readEnd1 + ", refEnd2=" + refEnd2 + ", cigar=" + cigarString() + ", flag=" + flag
				+ "]";
	}
}
