package au.edu.uq.imb.memesuite.db;

import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * A tgene database entry.
 */
public class TgeneDB implements JsonWr.JsonValue {
  private long listingId;
  private String name;
  private String description;
  private String genome_release;
  private String rna_source;
  private String tissues;
  private String histone_root;
  private String histones;
  private String max_link_distances;
  private String expression_root;
  private String annotation_file_name;
  private String transcript_types;
  private String use_gene_ids;
  private double lecat;

  public TgeneDB(
    long listingId,
    String name,
    String description,
    String genome_release,
    String rna_source,
    String tissues,
    String histone_root,
    String histones,
    String max_link_distances,
    String expression_root,
    String annotation_file_name,
    String transcript_types,
    String use_gene_ids,
    double lecat 
  ) {
    this.listingId = listingId;
    this.name = name;
    this.description = description;
    this.genome_release = genome_release;
    this.rna_source = rna_source;
    this.tissues = tissues;
    this.histone_root = histone_root;
    this.histones = histones;
    this.max_link_distances = max_link_distances;
    this.expression_root = expression_root;
    this.annotation_file_name = annotation_file_name;
    this.transcript_types = transcript_types;
    this.use_gene_ids = use_gene_ids;
    this.lecat = lecat;
  }

  public long getListingId() {
    return listingId;
  }
  public String getName() {
    return name;
  }
  public String getDescription() {
    return description;
  }
  public String getGenomeRelease() {
    return genome_release;
  }
  public String getRnaSource() {
    return rna_source;
  }
  public String getTissues() {
    return tissues;
  }
  public String getHistoneRoot() {
    return histone_root;
  }
  public String getHistones() {
    return histones;
  }
  public String getMaxLinkDistances() {
    return max_link_distances;
  }
  public String getExpressionRoot() {
    return expression_root;
  }
  public String getAnnotationFileName() {
    return annotation_file_name;
  }
  public String getTranscriptTypes() {
    return transcript_types;
  }
  public String getUseGeneIds() {
    return use_gene_ids;
  }
  public double getLecat() {
    return lecat;
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    out.property("name", getName());
    out.property("description", getDescription());
    out.property("genome_release", getGenomeRelease());
    out.property("tissues", getTissues());
    out.property("rna_source", getRnaSource());
    out.property("histones", getHistones());
    out.property("max_link_distances", getMaxLinkDistances());
    out.property("lecat", getLecat());
    out.property("use_gene_ids", getUseGeneIds());
    out.endObject();
  }
}
