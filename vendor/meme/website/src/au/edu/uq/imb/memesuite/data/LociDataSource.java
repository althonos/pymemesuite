package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.File;
import java.io.IOException;

/**
 * This class describes a loci file as a data source
 */
public class LociDataSource extends NamedFileDataSource implements LociInfo {
  private LociStats stats;

  /**
   * Create a loci data source from a file, the name of the file that the
   * submitter used and some basic information about the locis in the file.
   */
  public LociDataSource(File file, FileCoord.Name name,
      LociStats stats) {
    super(file, name);
    this.stats = stats;
  }

  @Override
  public long getLociCount() {
    return this.stats.getLociCount();
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    if (getOriginalName() != null) {
      out.property("source", "file");
      out.property("safe-file", getName());
      out.property("orig-file", getOriginalName());
    } else {
      out.property("source", "text");
    }
    out.property("count", getLociCount());
    out.endObject();
  }
}
