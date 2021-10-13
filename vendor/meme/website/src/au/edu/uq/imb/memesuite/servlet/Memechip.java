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


public class Memechip extends SubmitJob<Memechip.Data> {
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
  private ComponentAdvancedOptions memeOpts;
  private ComponentAdvancedOptions stremeOpts;
  private ComponentAdvancedOptions centrimoOpts;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.memechip");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public AlphabetDataSource alphabet;
    public SequenceDataSource sequences;
    public SequenceDataSource negSeq;
    public MotifInfo motifs;
    public String disc_mode;
    // universal options
    public boolean dna2rna;
    int minw;
    int maxw;
    public Background background;
    // MEME options
    String memeOptMode;
    int memeOptNMotifs;
    Integer memeOptMinSites;
    Integer memeOptMaxSites;
    boolean memeOptPal; 	// look for palindromes
    boolean memeOptNorand;	// MEME won't randomize sequence order
    // STREME options
    double stremeOptP;
    Integer stremeOptM;
    // CentriMo options
    double centrimoOptScore;
    Integer centrimoOptMaxReg;
    double centrimoOptEThresh;
    boolean centrimoOptLocal;
    boolean centrimoOptStoreIds;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.property("sequences", sequences);
      if (negSeq != null) out.property("negSeq", negSeq);
      out.property("disc_mode", disc_mode);
      out.property("motifs", motifs);
      out.property("dna2rna", dna2rna);
      out.property("minw", minw);
      out.property("maxw", maxw);
      out.property("background", background);
      out.property("memeOptMode", memeOptMode);
      out.property("memeOptNMotifs", memeOptNMotifs);
      if (memeOptMinSites != null) out.property("memeOptMinSites", memeOptMinSites);
      if (memeOptMaxSites != null) out.property("memeOptMaxSites", memeOptMaxSites);
      out.property("memeOptPal", memeOptPal);
      out.property("memeOptNorand", memeOptNorand);
      if (stremeOptP > 0) out.property("stremeOptP", stremeOptP);
      if (stremeOptM != null) out.property("stremeOptM", stremeOptM);
      out.property("centrimoOptScore", centrimoOptScore);
      if (centrimoOptMaxReg != null) out.property("centrimoOptMaxReg", centrimoOptMaxReg);
      out.property("centrimoOptEThresh", centrimoOptEThresh);
      out.property("centrimoOptLocal", centrimoOptLocal);
      out.property("centrimoOptStoreIds", centrimoOptStoreIds);
      out.endObject();
    }
    
    @Override
    public String email() {
      return email;
    }
  
    @Override
    public String description() {
      return description;  // generated code
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
        addArgs(args, "-alpha", sequences.guessAlphabet().name());
      }
      if (dna2rna) addArgs(args, "-dna2rna");
      if (motifs instanceof MotifDataSource) {
        addArgs(args, "-upmotif", ((MotifDataSource) motifs).getName());
      }
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bfile", background.getBfile().getName());
      } else {
        addArgs(args, "-order", background.getSource().getGeneratedOrder());
      }
      if (disc_mode == "psp") addArgs(args, "-psp-gen");
      if (negSeq != null) addArgs(args, "-neg", negSeq.getName());
      addArgs(args, "-minw", minw);
      addArgs(args, "-maxw", maxw);
      // MEME specific arguments
      addArgs(args, "-meme-mod", memeOptMode);
      addArgs(args, "-meme-nmotifs", memeOptNMotifs);
      if (memeOptMinSites != null) addArgs(args, "-meme-minsites", memeOptMinSites);
      if (memeOptMaxSites != null) addArgs(args, "-meme-maxsites", memeOptMaxSites);
      if (memeOptPal) addArgs(args, "-meme-pal");
      if (memeOptNorand) addArgs(args, "-meme-norand");
      // STREME specific options
      if (stremeOptP > 0) addArgs(args, "-streme-pvt", stremeOptP);
      if (stremeOptM != null) addArgs(args, "-streme-nmotifs", stremeOptM);
      // CentriMo specific options
      if (centrimoOptLocal) addArgs(args, "-centrimo-local");
      addArgs(args, "-centrimo-score", centrimoOptScore);
      if (centrimoOptMaxReg != null) addArgs(args, "-centrimo-maxreg", centrimoOptMaxReg);
      addArgs(args, "-centrimo-ethresh", centrimoOptEThresh);
      if (!centrimoOptStoreIds) addArgs(args, "-centrimo-noseq");
      // sequences
      addArgs(args, sequences.getName());
      // motif databases
      if (motifs instanceof MotifDB) {
        for (MotifDBFile dbFile : ((MotifDB) motifs).getMotifFiles()) {
          addArgs(args, dbFile.getFileName());
        }
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (alphabet != null) list.add(alphabet);
      list.add(sequences);
      if (negSeq != null) list.add(negSeq);
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
      if (negSeq != null) {
        if (!negSeq.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              negSeq.getFile());
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

  public Memechip() {
    super("MEMECHIP", "MEME-ChIP");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/meme-chip.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/meme-chip_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    control = new ComponentSequences(context, tmplMain.getSubtemplate("control"));
    background = new ComponentBfile(context, tmplMain.getSubtemplate("bfile"));
    jobDetails = new ComponentJobDetails(cache);
    universalOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("universal_opts"));
    memeOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("meme_opts"));
    stremeOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("streme_opts"));
    centrimoOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("centrimo_opts"));
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
    main.set("meme_opts", memeOpts.getComponent());
    main.set("streme_opts", stremeOpts.getComponent());
    main.set("centrimo_opts", centrimoOpts.getComponent());
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
    FileCoord.Name negSeqName = namer.createName("control_sequences.fa");
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
    // get the sequences
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
    // get the discovery mode
    data.disc_mode = paramChoice(request, "disc_mode", "classic", "de", "psp");
    // get any negative sequences sequences
    if (data.disc_mode != "classic") {
      data.negSeq = (SequenceDataSource)control.getSequences(alph, negSeqName, request, feedback);
      if (data.sequences.getCustomName() != null) data.negSeq.setCustomName(data.sequences.getCustomName());
    }
    // get the motifs
    data.motifs = motifs.getMotifs(alph, motifsName, request, feedback);
    // get the universal options
    data.background = background.getBfile(backgroundName, request, feedback);
    data.dna2rna = paramBool(request, "dna2rna");
    data.minw = paramInteger(feedback, "minimum motif width", request, "minw", 3, 30, 6);
    data.maxw = paramInteger(feedback, "maximum motif width", request, "maxw", 3, 30, 15);
    if (data.minw > data.maxw) {
      feedback.whine("The minimum motif width is larger than the maximum motif width.");
    }
    // get the MEME options
    data.memeOptMode = paramChoice(request, "meme_dist", "zoops", "oops", "anr");
    data.memeOptNMotifs = paramInteger(feedback, "MEME number of motifs", request, "meme_nmotifs", 0, null, 3);
    if (paramBool(request, "meme_minsites_enable")) {
      data.memeOptMinSites = paramInteger(feedback, "MEME minimum motif sites", request, "meme_minsites", 2, 600, 2);
    } else {
      data.memeOptMinSites = null;
    }
    if (paramBool(request, "meme_maxsites_enable")) {
      data.memeOptMaxSites = paramInteger(feedback, "MEME maximum motif sites", request, "meme_maxsites", 2, 600, 600);
      if (data.memeOptMinSites != null && data.memeOptMinSites > data.memeOptMaxSites) {
        feedback.whine("The MEME minimum motif sites is larger than the MEME maximum motif sites.");
      }
    } else {
      data.memeOptMaxSites = null;
    }
    data.memeOptPal = paramBool(request, "meme_pal");
    data.memeOptNorand = paramBool(request, "meme_norand");
    // get the STREME options
    if ("streme_enable_nmotifs".equals(paramRequire(request, "streme_srch_limit"))) {
      // motif count
      data.stremeOptM = paramInteger(feedback, "STREME motif count", request, "streme_nmotifs", 0, null, 5);
      data.stremeOptP = 0;
    } else {
      // p-value threshold
      data.stremeOptP = paramNumber(feedback, "STREME p-value threshold", request, "streme_pthresh", 0.0, null, 0.05);
      data.stremeOptM = null;
    }
    // get the CentriMo options
    data.centrimoOptScore = paramNumber(feedback, "CentriMo match score threshold", request, "centrimo_score", null, null, 5.0);
    if (paramBool(request, "centrimo_maxreg_enable")) {
      data.centrimoOptMaxReg = paramInteger(feedback, "CentriMo maximum region", request, "centrimo_maxreg", 0, null, 200);
    } else {
      data.centrimoOptMaxReg = null;
    }
    data.centrimoOptEThresh = paramNumber(feedback, "CentriMo E-value threshold", request, "centrimo_ethresh", 0.0, null, 10.0);
    data.centrimoOptLocal = paramBool(request, "centrimo_local");
    data.centrimoOptStoreIds = paramBool(request, "centrimo_store_ids");
    return data;
  }

}

