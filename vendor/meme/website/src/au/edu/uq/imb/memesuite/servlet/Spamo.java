package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.MotifDB;
import au.edu.uq.imb.memesuite.db.MotifDBFile;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
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

public class Spamo extends SubmitJob<Spamo.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentMotifs primary;
  private ComponentMotifs secondaries;
  private ComponentSequences sequences;
  private ComponentBfile background;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advBtn;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.spamo");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public SequenceDataSource sequences;
    public MotifDataSource primary;
    public MotifInfo secondaries;
    public boolean dumpseqs;
    public boolean xalph;
    public Integer margin;
    public Background background;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.property("sequences", sequences);
      out.property("primary", primary);
      out.property("secondaries", secondaries);
      out.property("dumpseqs", dumpseqs);
      out.property("xalph", xalph);
      out.property("background", background);
      if (margin != null) out.property("margin", margin);
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
  //    spamo_webservice [options] <sequences file> <primary motif> <secondary db patterns>
  //
  //      Options:
  //        -uploaded <file>  file containing uploaded secondary motif database
  //        -margin <margin>  margin parameter passed to spamo
  //        -bgfile <file>    background file (0-order)
  //        -dumpseqs         dump the sequence matches to a file for each significant primary/secondary
  //        -help             brief help message
      StringBuilder args = new StringBuilder();
      if (secondaries instanceof MotifDataSource) {
        addArgs(args, "-uploaded", ((MotifDataSource) secondaries).getName());
      }
      if (margin != null) addArgs(args, "-margin", margin);
      if (dumpseqs) addArgs(args, "-dumpseqs");
      if (xalph) addArgs(args, "-xalph");
      addArgs(args, sequences.getName());
      addArgs(args, primary.getName());
      if (secondaries instanceof MotifDB) {
        for (MotifDBFile file : ((MotifDB) secondaries).getMotifFiles()) {
          addArgs(args, file.getFileName());
        }
      }
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bgfile", background.getBfile().getName());
      } else if (background.getSource() == Background.Source.UNIFORM) {
        addArgs(args, "-bgfile", "--uniform--"); // Use uniform background
      } else if (background.getSource() == Background.Source.MEME) {
        addArgs(args, "-bgfile", "--motif--"); // Use background in motif file
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (sequences != null) list.add(sequences);
      if (primary != null) list.add(primary);
      if (secondaries != null && secondaries instanceof MotifDataSource) {
        list.add((MotifDataSource) secondaries);
      }
      if (background.getSource() == Background.Source.FILE)
        list.add(background.getBfile());
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (sequences != null) {
        if (!sequences.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              sequences.getFile());
        }
      }
      if (primary != null) {
        if (!primary.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              primary.getFile());
        }
      }
      if (secondaries != null && secondaries instanceof MotifDataSource) {
        if (!((MotifDataSource) secondaries).getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              ((MotifDataSource) secondaries).getFile());
        }
      }
      if (background.getSource() == Background.Source.FILE) {
        if (!background.getBfile().getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              background.getBfile().getFile());
        }
      }
    }
  }

  public Spamo() {
    super("SPAMO", "SpaMo");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/spamo.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/spamo_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    primary = new ComponentMotifs(context, tmplMain.getSubtemplate("primary"));
    secondaries = new ComponentMotifs(context, tmplMain.getSubtemplate("secondaries"));
    background = new ComponentBfile(context, tmplMain.getSubtemplate("bfile"));
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
    main.set("help", new HTMLSub[]{header.getHelp(), primary.getHelp(),
        sequences.getHelp(), jobDetails.getHelp(), advBtn.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("sequences", sequences.getComponent());
    main.set("primary", primary.getComponent(request.getParameter("motifs_embed")));
    main.set("secondaries", secondaries.getComponent());
    main.set("bfile", background.getComponent());
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
    FileCoord.Name sequencesName = namer.createName("sequences.fa");
    FileCoord.Name primaryName = namer.createName("primary.meme");
    FileCoord.Name secondariesName = namer.createName("secondaries.meme");
    FileCoord.Name backgroundName = namer.createName("background");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    Alph secondaryAlph = null;
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the motifs
    data.primary =  (MotifDataSource)primary.getMotifs(primaryName, request, feedback);
    if (data.primary != null) alph = data.primary.getAlphabet();
    data.secondaries = secondaries.getMotifs(secondariesName, request, feedback);
    if (data.secondaries != null) secondaryAlph = data.secondaries.getAlphabet();
    if (alph != null && secondaryAlph != null && !alph.equals(secondaryAlph)) {
      if (Alph.checkCoreSubset(secondaryAlph, alph) == 0) {
        feedback.whine("The target motifs have the " + secondaryAlph +
            " alphabet which is not a subset of the " + alph + " alphabet.");
      }
      data.xalph = true;
    } else {
      data.xalph = false;
    }
    // get the sequences
    data.sequences = (SequenceDataSource)sequences.getSequences(alph, sequencesName, request, feedback);
    // get the options
    data.background = background.getBfile(backgroundName, request, feedback);
    data.dumpseqs = paramBool(request, "dumpseqs");
    data.margin = null; //TODO
    return data;
  }

}

