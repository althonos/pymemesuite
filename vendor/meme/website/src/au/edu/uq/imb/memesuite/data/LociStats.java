package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.JsonWr;
import au.edu.uq.imb.memesuite.util.SampleStats;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

/**
 * Calculate statistics on the loci
 */
public final class LociStats implements LociInfo {
  private boolean frozen;
  private SampleStats stats;

  /**
   * Class constructor
   */
  public LociStats() {
    this.stats = new SampleStats(true);
  }

  /**
   * Return the total count of loci
   * @return The total count of loci
   */
  public long getLociCount() {
    //return this.stats.getCount();
    return 0;
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    out.property("type", "loci");
    out.property("count", getLociCount());
    out.endObject();
  }
}
