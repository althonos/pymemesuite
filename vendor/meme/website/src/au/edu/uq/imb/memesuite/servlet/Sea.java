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

public class Sea extends SubmitJob<Sea.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentSequenceAlphabet alphabet;
  private ComponentSequences sequences;
  private ComponentSequences control;
  private ComponentMotifs motifs;
  private ComponentBfile background;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advancedOptions;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public AlphabetDataSource alphabet;
    public SequenceDataSource posSeq;
    public SequenceDataSource negSeq;
    public MotifInfo motifs;
    public Background background;
    public double evt;
    public Integer order;
    public String align;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      if (alphabet != null) out.property("alphabet", alphabet);
      out.property("posSeq", posSeq);
      if (negSeq != null) out.property("negSeq", negSeq);
      out.property("motifs", motifs);
      out.property("background", background);
      if (evt > 0) out.property("evt", evt);
      if (order != null) out.property("order", order);
      out.property("align", align);
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
      addArgs(args, "-p", posSeq.getName());
      if (negSeq != null) addArgs(args, "-n", negSeq.getName());
      if (alphabet != null) { addArgs(args, "-xalph", alphabet.getName()); }
      if (evt > 0) addArgs(args, "-thresh", evt);
      if (order != null) addArgs(args, "-order", order);
      if (align != null) addArgs(args, "-align", align);
      if (motifs instanceof MotifDataSource) {
        addArgs(args, "-m", ((MotifDataSource) motifs).getName());
      } else if (motifs instanceof MotifDB) {
        for (MotifDBFile file : ((MotifDB) motifs).getMotifFiles()) {
          addArgs(args, "-m", "db/" + file.getFileName());
        }
      }
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bfile", background.getBfile().getName());	// use bfile
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> sources = new ArrayList<DataSource>();
      if (alphabet != null) sources.add(alphabet);
      if (posSeq != null) sources.add(posSeq);
      if (negSeq != null) sources.add(negSeq);
      if (motifs instanceof MotifDataSource) sources.add((MotifDataSource) motifs);
      if (background.getSource() == Background.Source.FILE) sources.add(background.getBfile());
      return sources;
    }
  
    @Override
    public void cleanUp() {
      if (alphabet != null) alphabet.getFile().delete();
      if (posSeq != null) posSeq.getFile().delete();
      if (negSeq != null) negSeq.getFile().delete();
      if (motifs != null && motifs instanceof MotifDataSource) ((MotifDataSource) motifs).getFile().delete();
      if (background != null && background.getSource() == Background.Source.FILE) background.getBfile().getFile().delete();
    }
  }

  public Sea() {
    super("SEA", "SEA");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the template
    this.tmplMain = cache.loadAndCache("/WEB-INF/templates/sea.tmpl");
    this.tmplVerify = cache.loadAndCache("/WEB-INF/templates/sea_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    control = new ComponentSequences(context, tmplMain.getSubtemplate("control"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    background = new ComponentBfile(context, tmplMain.getSubtemplate("bfile"));
    jobDetails = new ComponentJobDetails(cache);
    advancedOptions = new ComponentAdvancedOptions(cache);
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
    main.set("help", new HTMLSub[]{header.getHelp(), alphabet.getHelp(), sequences.getHelp(),
        motifs.getHelp(), background.getHelp(), jobDetails.getHelp(), advancedOptions.getHelp(), 
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("alphabet", alphabet.getComponent());
    main.set("sequences", sequences.getComponent());
    main.set("control", control.getComponent());
    main.set("bfile", background.getComponent());
    main.set("motifs", motifs.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("advanced_options", advancedOptions.getComponent());
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
      // get the positive sequences
      data.posSeq = (SequenceDataSource)sequences.getSequences(alph, posSeqName, request, feedback);
      if (data.posSeq != null) {
        if (alph == null) {
          // Get the alphabet from the positive sequences by guessing.
          alph = data.posSeq.guessAlphabet().getAlph();
        } else {
          // Custom alphabet given; set its name in the posSeq object.
          data.posSeq.setCustomName(alph.getName());
        }
      }
      // get the discriminative sequences
      if (paramBool(request, "discr")) {
        data.negSeq = (SequenceDataSource)control.getSequences(alph, negSeqName, request, feedback);
        if (data.posSeq.getCustomName() != null) data.negSeq.setCustomName(data.posSeq.getCustomName());
      }
      // get the motifs
      data.motifs = motifs.getMotifs(alph, motifsName, request, feedback);
      // E-value threshold
      data.evt = paramNumber(feedback, "<i>E</i>-value threshold", request,
          "evt", 0.0, 1e300, 10.0);
      // Order
      if (paramBool(request, "order_enable")) {
        data.order = paramInteger(feedback, "Markov order", request, "order", 0, 4, 0);
      }
      // Background model
      data.background = background.getBfile(backgroundName, request, feedback);
      error = false;
      // Align
      data.align = paramChoice(request, "align", "left", "center", "right");
    } finally {
      if (error) data.cleanUp();
    }
    return data;
  }
}
