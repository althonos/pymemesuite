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

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;

public class Xstreme extends SubmitJob<Xstreme.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentMotifs motifs;
  private ComponentSequenceAlphabet alphabet;
  private ComponentSequences sequences;
  private ComponentSequences control;
  private ComponentBfile background;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions universalOpts;
  private ComponentAdvancedOptions stremeOpts;
  private ComponentAdvancedOptions memeOpts;
  private ComponentAdvancedOptions seaOpts;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public AlphabetDataSource alphabet;
    public SequenceDataSource posSeq;
    public SequenceDataSource negSeq;
    public boolean dna2rna;
    public MotifInfo motifs;
    // universal options
    public double evt;
    public Integer minw;
    public Integer maxw;
    public Background background;
    public Integer order;
    public Integer ctrim;
    public String align;
    // STREME options
    double stremeOptE;
    Integer stremeOptM;
    // MEME options
    double memeOptE;
    Integer memeOptM;
    String memeOptMode;
    // SEA options
    public boolean seaSeqs;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      if (alphabet != null) out.property("alphabet", alphabet);
      out.property("posSeq", posSeq);
      if (negSeq != null) out.property("negSeq", negSeq);
      out.property("dna2rna", dna2rna);
      out.property("motifs", motifs);
      if (evt > 0) out.property("evt", evt);
      out.property("minw", minw);
      out.property("maxw", maxw);
      if (order != null) out.property("order", order);
      if (ctrim != null) out.property("ctrim", ctrim);
      if (align != null) out.property("align", align);
      out.property("background", background);
      if (stremeOptE > 0) out.property("stremeOptE", stremeOptE);
      if (stremeOptM != null) out.property("stremeOptM", stremeOptM);
      if (memeOptE > 0) out.property("memeOptE", memeOptE);
      if (memeOptM != null) out.property("memeOptM", memeOptM);
      if (memeOptM != null) out.property("memeOptM", memeOptM);
      if (memeOptMode != null) out.property("memeOptMode", memeOptMode);
      out.property("seaSeqs", seaSeqs);
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
        addArgs(args, "-alphf", alphabet.getName());
      } else {
        String guess_alph = posSeq.guessAlphabet().name();
        addArgs(args, "-alpha", guess_alph);
      }
      addArgs(args, "-p", posSeq.getName());
      if (negSeq != null) addArgs(args, "-n", negSeq.getName());
      if (dna2rna) addArgs(args, "-dna2rna");
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bfile", background.getBfile().getName());
      }
      if (evt > 0) addArgs(args, "-evt", evt);
      addArgs(args, "-minw", minw);
      addArgs(args, "-maxw", maxw);
      // STREME specific options
      if (stremeOptE > 0) addArgs(args, "-streme-evt", stremeOptE);
      if (stremeOptM != null) addArgs(args, "-streme-nmotifs", stremeOptM);
      // MEME specific arguments
      if (memeOptE > 0) addArgs(args, "-meme-evt", memeOptE);
      if (memeOptM != null) addArgs(args, "-meme-nmotifs", memeOptM);
      if (memeOptMode != null) addArgs(args, "-meme-mod", memeOptMode);
      if (order != null) addArgs(args, "-order", order);
      if (ctrim != null) addArgs(args, "-ctrim", ctrim);
      if (align != null) addArgs(args, "-align", align);
      if (! seaSeqs) addArgs(args, "-sea-noseqs");
      // uploaded motifs
      if (motifs instanceof MotifDataSource) {
        addArgs(args, "-upmotif", ((MotifDataSource) motifs).getName());
      }
      // motif databases (no option switch)
      if (motifs instanceof MotifDB) {
        for (MotifDBFile dbFile : ((MotifDB) motifs).getMotifFiles()) {
          addArgs(args, dbFile.getFileName());
        }
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      List<DataSource> sources = new ArrayList<DataSource>();
      if (alphabet != null) sources.add(alphabet);
      if (posSeq != null) sources.add(posSeq);
      if (negSeq != null) sources.add(negSeq);
      if (motifs instanceof MotifDataSource) {
        sources.add((MotifDataSource) motifs);
      }
      if (background.getSource() == Background.Source.FILE) {
        sources.add(background.getBfile());
      }
      return sources;
    }
  
    @Override
    public void cleanUp() {
      if (alphabet != null) alphabet.getFile().delete();
      if (posSeq != null) posSeq.getFile().delete();
      if (negSeq != null) negSeq.getFile().delete();
      if (motifs != null & motifs instanceof MotifDataSource) ((MotifDataSource) motifs).getFile().delete();
      if (background.getSource() == Background.Source.FILE) background.getBfile().getFile().delete();
    }
  }

  public Xstreme() {
    super("XSTREME", "XSTREME");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    this.tmplMain = cache.loadAndCache("/WEB-INF/templates/xstreme.tmpl");
    this.tmplVerify = cache.loadAndCache("/WEB-INF/templates/xstreme_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    control = new ComponentSequences(context, tmplMain.getSubtemplate("control"));
    background = new ComponentBfile(context, tmplMain.getSubtemplate("bfile"));
    jobDetails = new ComponentJobDetails(cache);
    universalOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("universal_opts"));
    stremeOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("streme_opts"));
    memeOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("meme_opts"));
    seaOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("sea_opts"));
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
    HTMLSub main = this.tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(), alphabet.getHelp(), motifs.getHelp(),
        sequences.getHelp(), background.getHelp(), jobDetails.getHelp(), 
        universalOpts.getHelp(), submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("alphabet", alphabet.getComponent());
    main.set("motifs", motifs.getComponent());
    main.set("sequences", sequences.getComponent());
    main.set("control", control.getComponent());
    main.set("bfile", background.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("universal_opts", universalOpts.getComponent());
    main.set("streme_opts", stremeOpts.getComponent());
    main.set("meme_opts", memeOpts.getComponent());
    main.set("sea_opts", seaOpts.getComponent());
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
    FileCoord.Name posSeqName = namer.createName("sequences.fa");
    FileCoord.Name negSeqName = namer.createName("control_sequences.fa");
    FileCoord.Name motifsName = namer.createName("motifs.meme");
    FileCoord.Name backgroundName = namer.createName("background");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    // create the job data
    Data data =  new Data();
    boolean error = true;
    try {
      // get the email
      data.email = jobDetails.getEmail(request, feedback);
      // get the description
      data.description = jobDetails.getDescription(request);
      // get the alphabet
      data.alphabet = alphabet.getAlphabet(alphName, request, feedback);
      if (data.alphabet != null) alph = data.alphabet.getAlph();
      // get the primary sequences
      data.posSeq = (SequenceDataSource)sequences.getSequences(alph, posSeqName, request, feedback);
      if (data.posSeq != null) {
        if (alph == null) {
          // Get the alphabet from the primary sequences by guessing.
          alph = data.posSeq.guessAlphabet().getAlph();
        } else {
          // Custom alphabet given; set its name in the posSeq object.
          data.posSeq.setCustomName(alph.getName());
        }
      }
      // get the control sequences
      if (paramBool(request, "discr")) {
        data.negSeq = (SequenceDataSource)control.getSequences(alph, negSeqName, request, feedback);
        if (data.posSeq.getCustomName() != null) data.negSeq.setCustomName(data.posSeq.getCustomName());
      }
      data.dna2rna = paramBool(request, "dna2rna");
      // get the motifs
      data.motifs = motifs.getMotifs(alph, motifsName, request, feedback);
      // get the universal options
      data.evt = paramNumber(feedback, "E-value threshold", request, "evt", 0.0, 10.0, 0.05);
      data.minw = paramInteger(feedback, "Minimum motif width", request, "minw", 3, 30, 8);
      data.maxw = paramInteger(feedback, "Maximum motif width", request, "maxw", 3, 30, 15);
      if (data.minw > data.maxw) {
	feedback.whine("The minimum motif width is larger than the maximum motif width.");
      }
      if (paramBool(request, "order_enable")) {
        data.order = paramInteger(feedback, "Markov order", request, "order", 0, 4, 0);
      }
      if (paramBool(request, "ctrim_enable")) {
        data.ctrim = paramInteger(feedback, "Central trimming size", request, "ctrim", 2, null, 100);
      }
      // Align
      data.align = paramChoice(request, "align", "left", "center", "right");
      data.background = background.getBfile(backgroundName, request, feedback);
      // get the STREME options
      if ("streme_enable_ethresh".equals(paramOption(request, "streme_srch_limit", "default"))) {
	data.stremeOptE = paramNumber(feedback, "STREME E-value threshold", request, "streme_ethresh", 0.0, null, 0.05);
	data.stremeOptM = null;
      } else if ("streme_enable_nmotifs".equals(paramOption(request, "streme_srch_limit", "default"))) {
	data.stremeOptM = paramInteger(feedback, "STREME number of motifs", request, "streme_nmotifs", 0, null, 5);
	data.stremeOptE = 0;
      }
      // get the MEME options
      if ("meme_enable_ethresh".equals(paramOption(request, "meme_srch_limit", "default"))) {
	data.memeOptE = paramNumber(feedback, "MEME E-value threshold", request, "meme_ethresh", 0.0, null, 0.05);
	data.memeOptM = null;
      } else if ("meme_enable_nmotifs".equals(paramOption(request, "meme_srch_limit", "default"))) {
	data.memeOptM = paramInteger(feedback, "MEME number of motifs", request, "meme_nmotifs", 0, null, 3);
	data.memeOptE = 0;
      }
      data.memeOptMode = paramChoice(request, "meme_dist", "zoops", "oops", "anr");
      // get the SEA options
      data.seaSeqs = paramBool(request, "sea_seqs");
      error = false;
    } finally {
      if (error) data.cleanUp();
    }
    return data;
  }
}
