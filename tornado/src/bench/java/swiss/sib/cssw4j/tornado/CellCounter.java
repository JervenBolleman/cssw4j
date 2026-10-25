package swiss.sib.cssw4j.tornado;

import org.openjdk.jmh.annotations.AuxCounters;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;

/**
 * Counts dynamic-programming cells for benchmarks whose work per call depends on a parameter. JMH reports it in
 * throughput mode as a secondary result {@code :cells} in ops per time unit; cells per ms divided by 10^6 is GCUPS.
 */
@State(Scope.Thread)
@AuxCounters(AuxCounters.Type.OPERATIONS)
public class CellCounter {
	public long cells;

	@Setup(Level.Iteration)
	public void reset() {
		cells = 0;
	}
}
