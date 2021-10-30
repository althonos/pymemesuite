package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.File;
import java.io.IOException;

/**
 * This class describes a sequence file as a data source
 */
public class SequenceDataSource extends NamedFileDataSource implements SequenceInfo {
  private SequenceStats stats;
  private String customName;

  /**
   * Create a sequence data source from a file, the name of the file that the
   * submitter used and some basic information about the sequences in the file.
   */
  public SequenceDataSource(File file, FileCoord.Name name,
      SequenceStats stats) {
    super(file, name);
    this.stats = stats;
    this.customName = null;
  }

  /* Set the alphabet name because it is a custom alphabet name */
  public void setCustomName(String customName) {
    this.customName = customName;
  }

  /* Get the custom alphabet */
  public String getCustomName() {
    return this.customName;
  }

  @Override
  public AlphStd guessAlphabet() {
    return this.stats.guessAlphabet();
  }

  @Override
  public boolean checkAlphabet(Alph alph) {
    return this.stats.checkAlphabet(alph);
  }

  @Override
  public long getSequenceCount() {
    return this.stats.getSequenceCount();
  }

  @Override
  public long getTotalLength() {
    return this.stats.getTotalLength();
  }

  @Override
  public long getMinLength() {
    return this.stats.getMinLength();
  }

  @Override
  public long getMaxLength() {
    return this.stats.getMaxLength();
  }

  @Override
  public double getAverageLength() {
    return this.stats.getAverageLength();
  }

  @Override
  public double getStandardDeviationLength() {
    return this.stats.getStandardDeviationLength();
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
    out.property("alphabet", guessAlphabet().getAlph());
    out.property("custom_name", this.customName);
    out.property("count", getSequenceCount());
    out.property("min", getMinLength());
    out.property("max", getMaxLength());
    out.property("avg", getAverageLength());
    out.property("total", getTotalLength());
    out.endObject();
  }
}
