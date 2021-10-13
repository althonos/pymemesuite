package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.LociDataSource;
import au.edu.uq.imb.memesuite.db.TgeneDB;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;


public class Tgene extends SubmitJob<Tgene.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentLoci loci;
  private ComponentTgene tgenePanel;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advBtn;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.tgene");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public LociDataSource loci;
    public TgeneDB tgenePanel;
    public double max_pvalue;
    public boolean closest_locus;
    public boolean closest_tss;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      boolean have_panel = ! tgenePanel.getTissues().equals("none");
      out.startObject();
      out.property("loci", loci);
      out.property("annotation_file_name", tgenePanel.getAnnotationFileName());
      out.property("transcript_types", tgenePanel.getTranscriptTypes().replaceAll(",", ", "));
      out.property("max_link_distances", tgenePanel.getMaxLinkDistances().replaceAll(",", ", "));
      out.property("max_pvalue", max_pvalue);
      if (have_panel) {
        out.property("panel_name", tgenePanel.getName());
        out.property("panel_description", tgenePanel.getDescription());
        out.property("genome_release", tgenePanel.getGenomeRelease());
	out.property("tissues", tgenePanel.getTissues().replaceAll(",", ", "));
	out.property("histones", tgenePanel.getHistones().replaceAll(",", ", "));
	out.property("rna_source", tgenePanel.getRnaSource());
	out.property("expression_root", tgenePanel.getExpressionRoot());
	out.property("use_gene_ids", tgenePanel.getUseGeneIds().equals("true"));
	out.property("lecat", tgenePanel.getLecat());
      } else {
        out.property("genome_name", tgenePanel.getName());
        out.property("genome_description", tgenePanel.getDescription());
        out.property("genome_release", tgenePanel.getGenomeRelease());
      }
      out.property("closest_locus", closest_locus);
      out.property("closest_tss", closest_tss);
      out.endObject();
    }

    @Override
    public String email() {
      return email;
    }
  
    @Override
    public String description() {
      return description;
    }

    @Override
    public boolean immediateRun() {
      return false;
    }

    @Override
    public String emailTemplate() {
      return tmplVerify.getSubtemplate("message").toString();
    }

    @Override
    public String cmd() {
    // Refer to scripts/tgene_webservice.pl.in
      StringBuilder args = new StringBuilder();
      addArgs(args, loci.getName());
      addArgs(args, "-annotation-file-name", "db/"+tgenePanel.getAnnotationFileName());
      addArgs(args, "-transcript-types", tgenePanel.getTranscriptTypes());
      addArgs(args, "-max-link-distances", tgenePanel.getMaxLinkDistances());
      addArgs(args, "-max-pvalue", max_pvalue);
      if (tgenePanel.getTissues() != ".") {
	addArgs(args, "-tissues", tgenePanel.getTissues());
	addArgs(args, "-histone-root", "db/"+tgenePanel.getHistoneRoot());
	addArgs(args, "-histones", tgenePanel.getHistones());
	addArgs(args, "-rna-source", tgenePanel.getRnaSource());
	addArgs(args, "-expression-root", "db/"+tgenePanel.getExpressionRoot());
	if (tgenePanel.getUseGeneIds().equals("true")) { 
	  addArgs(args, "-use-gene-ids");
	};
	addArgs(args, "-lecat", tgenePanel.getLecat());
      }
      if (! closest_locus) addArgs(args, "-no-closest-locus");
      if (! closest_tss) addArgs(args, "-no-closest-tss");
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (loci != null) list.add(loci);
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (loci != null) {
        if (!loci.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              loci.getFile());
        }
      }
    }
  }

  public Tgene() {
    super("TGENE", "T-Gene");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplMain = cache.loadAndCache("/WEB-INF/templates/tgene.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/tgene_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    loci = new ComponentLoci(context, tmplMain.getSubtemplate("loci"));
    tgenePanel = new ComponentTgene(context, tmplMain.getSubtemplate("tgene_panel"));
    jobDetails = new ComponentJobDetails(cache);
    advBtn = new ComponentAdvancedOptions(cache);
    submitReset = new ComponentSubmitReset(cache, jobTable.getCount(), jobTable.getDuration());
    footer = new ComponentFooter(cache, msp);
  }

  @Override
  public String title() {
    return tmplVerify.getSubtemplate("title").toString();
  }

  @Override
  public String subtitle() {
    return tmplVerify.getSubtemplate("subtitle").toString();
  }

  @Override
  public String logoPath() {
    return tmplVerify.getSubtemplate("logo").toString();
  }

  @Override
  public String logoAltText() {
    return tmplVerify.getSubtemplate("alt").toString();
  }

  @Override
  protected void displayForm(HttpServletRequest request, HttpServletResponse response, long quotaMinWait) throws IOException {
    HTMLSub main = tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(), loci.getHelp(),
        tgenePanel.getHelp(), jobDetails.getHelp(), advBtn.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("loci", loci.getComponent());
    main.set("tgene_panel", tgenePanel.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("advanced_options", advBtn.getComponent());
    main.set("submit_reset", submitReset.getComponent(quotaMinWait));
    main.set("footer", footer.getComponent());
    response.setContentType("text/html; charset=UTF-8");
    main.output(response.getWriter());
  }

  @Override
  protected Data checkParameters(FeedbackHandler feedback,
      HttpServletRequest request) throws IOException, ServletException {
    // setup default file names
    FileCoord namer = new FileCoord();
    FileCoord.Name lociName = namer.createName("loci.bed");
    namer.createName("description");
    namer.createName("uuid");
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the loci
    data.loci = (LociDataSource)loci.getLoci(lociName, request, feedback);
    // get the panel 
    data.tgenePanel = tgenePanel.getTgenePanel(request);
    // get advanced options
    data.max_pvalue = paramNumber(feedback, "maximum p-value", request, "max_pvalue", 0.0, 1.0, 0.05);
    data.closest_locus = paramBool(request, "closest_locus");
    data.closest_tss = paramBool(request, "closest_tss");

    return data;
  } // checkParameters
}
