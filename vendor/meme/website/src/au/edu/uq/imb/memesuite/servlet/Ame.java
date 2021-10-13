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

public class Ame extends SubmitJob<Ame.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentSequenceAlphabet alphabet;
  private ComponentSequences sequences;
  private ComponentSequences control;
  private ComponentBfile background;
  private ComponentMotifs motifs;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advBtn;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.ame");
  
  protected class Data extends SubmitJob.JobData {
    public String control_type;
    public String email;
    public String description;
    public String method;
    public String scoring;
    public Double hit_lo_fraction;
    public Integer kmer;
    public Double evalue_report_threshold;
    public AlphabetDataSource alphabet;
    public SequenceDataSource sequences;
    public SequenceDataSource control;
    public Background background;
    public MotifInfo motifs;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      if (alphabet != null) out.property("alphabet", alphabet);
      out.property("sequences", sequences);
      if (control != null) out.property("control", control);
      out.property("motifs", motifs);
      out.property("background", background);
      out.property("method", method);
      out.property("scoring", scoring);
      if (control_type != null) out.property("control_type", control_type);
      if (hit_lo_fraction != null) out.property("hit_lo_fraction", hit_lo_fraction);
      if (kmer != null) out.property("kmer", kmer);
      if (evalue_report_threshold != null) out.property("evalue_report_threshold", evalue_report_threshold);
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
      StringBuilder args = new StringBuilder();
      if (alphabet != null) {
        addArgs(args, "-xalph", alphabet.getName());
      }
      addArgs(args, "-control_type", control_type);
      if (control != null) {
        addArgs(args, "-control", control.getName());
      }
      addArgs(args, "-method", method);
      addArgs(args, "-scoring", scoring);
      if (hit_lo_fraction != null) {
        addArgs(args, "-hit-lo-fraction", hit_lo_fraction);
      }
      if (kmer != null) {
        addArgs(args, "-kmer", kmer);
      }
      if (evalue_report_threshold != null) {
        addArgs(args, "-evalue-report-threshold", evalue_report_threshold);
      }
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bgfile", background.getBfile().getName());
      } else if (background.getSource() == Background.Source.UNIFORM) {
        addArgs(args, "-bgfile", "--uniform--"); // Use uniform background
      } else if (background.getSource() == Background.Source.MEME) {
        addArgs(args, "-bgfile", "--motif--"); // Use background in motif file
      }
      addArgs(args, sequences.getName());
      if (motifs instanceof MotifDataSource) {
        addArgs(args, ((MotifDataSource) motifs).getName());
      } else if (motifs instanceof MotifDB){
        for (MotifDBFile file : ((MotifDB) motifs).getMotifFiles()) {
          addArgs(args, "db/" + file.getFileName());
        }
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (alphabet != null) list.add(alphabet);
      list.add(sequences);
      if (control != null) list.add(control);
      if (motifs instanceof MotifDataSource) {
        list.add((MotifDataSource) motifs);
      }
      if (background.getSource() == Background.Source.FILE) {
        list.add(background.getBfile());
      }
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (alphabet != null) {
        if (!alphabet.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              alphabet.getFile());
        }
      }
      if (sequences != null) {
        if (!sequences.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              sequences.getFile());
        }
      }
      if (control != null) {
        if (!control.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              control.getFile());
        }
      }
      if (motifs != null && motifs instanceof MotifDataSource) {
        if (!((MotifDataSource) motifs).getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              ((MotifDataSource) motifs).getFile());
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

  public Ame() {
    super("AME", "AME");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/ame.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/ame_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    control = new ComponentSequences(context, tmplMain.getSubtemplate("control"));
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
    main.set("help", new HTMLSub[]{header.getHelp(), alphabet.getHelp(), motifs.getHelp(),
        sequences.getHelp(), jobDetails.getHelp(), advBtn.getHelp(), background.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("alphabet", alphabet.getComponent());
    main.set("motifs", motifs.getComponent());
    main.set("sequences", sequences.getComponent());
    main.set("control", control.getComponent());
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
    FileCoord.Name alphName = namer.createName("alphabet.alph");
    FileCoord.Name sequencesName = namer.createName("sequences.fa");
    FileCoord.Name controlName = namer.createName("control.fa");
    FileCoord.Name motifsName = namer.createName("motifs.meme");
    FileCoord.Name backgroundName = namer.createName("background");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the alphabet
    data.alphabet = alphabet.getAlphabet(alphName, request, feedback);
    if (data.alphabet != null) alph = data.alphabet.getAlph();
    // get the motifs
    data.motifs = motifs.getMotifs(alph, motifsName, request, feedback);
    // get the primary sequences
    data.sequences = (SequenceDataSource)sequences.getSequences(alph, sequencesName, request, feedback);
    if (data.sequences != null) {
      if (alph == null) {
        // Get the alphabet from the positive sequences by guessing.
        alph = data.sequences.guessAlphabet().getAlph();
      } else {
        // Custom alphabet given; set its name in the posSeq object.
        data.sequences.setCustomName(alph.getName());
      }
    }
    // get the control sequences
    data.control_type = WebUtils.paramChoice(request, "control_type", "shuffle", "user", "none");
    if (data.control_type == "user") {
      data.control = (SequenceDataSource)control.getSequences(alph, controlName, request, feedback);
      if (data.sequences.getCustomName() != null) data.control.setCustomName(data.sequences.getCustomName());
    } else {
      data.control = null;
    }
    // get the value of kmer
    if (data.control_type == "shuffle") {
      data.kmer = WebUtils.paramInteger(feedback, "shuffling preserves frequencies of words this size", request, "kmer", 1, 10, 2);
    } else {
      data.kmer = null;
   }
    data.method = WebUtils.paramChoice(request, "method", "fisher", "ranksum", "pearson", "spearman", "3dmhg", "4dmhg");
    data.scoring = WebUtils.paramChoice(request, "scoring", "avg", "max", "sum", "totalhits");
    data.hit_lo_fraction = WebUtils.paramNumber(feedback, "motif match log-odds fraction", request, "hit_lo_fraction", 0.0, 1.0, 0.25);

    // get advanced options
    data.evalue_report_threshold = paramNumber(feedback, "E-value report threshold", request, "evalue_report_threshold", 0.0, 1e300, 10.0);
    data.background = background.getBfile(backgroundName, request, feedback);

    return data;
  }
}
