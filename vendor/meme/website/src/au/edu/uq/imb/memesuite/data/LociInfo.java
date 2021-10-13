package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.JsonWr;

/**
 * Information about a group of loci.
 */
public interface LociInfo extends JsonWr.JsonValue {

  /**
   * Return the number of loci in the source.
   * @return The number of loci
   */
  public long getLociCount();

}
