package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.Alph;
import au.edu.uq.imb.memesuite.data.AlphabetDataSource;
import au.edu.uq.imb.memesuite.data.SequenceDataSource;
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


public class Dreme extends SubmitJob<Dreme.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentSequenceAlphabet alphabet;
  private ComponentSequences sequences;
  private ComponentSequences control;
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
    public double ethresh;
    public Integer nMotifs;
    public boolean norc;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      if (alphabet != null) out.property("alphabet", alphabet);
      out.property("posSeq", posSeq);
      if (negSeq != null) out.property("negSeq", negSeq);
      if (ethresh > 0) out.property("ethresh", ethresh);
      if (nMotifs != null) out.property("nMotifs", nMotifs);
      out.property("norc", norc);
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
        addArgs(args, "-alpha", posSeq.guessAlphabet().name());
      }
      if (norc) addArgs(args, "-norc");
      if (ethresh > 0) addArgs(args, "-e", ethresh);
      if (nMotifs != null) addArgs(args, "-m", nMotifs);
      if (negSeq != null) addArgs(args, "-n", negSeq.getName());
      addArgs(args, posSeq.getName());
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      List<DataSource> sources = new ArrayList<DataSource>();
      if (alphabet != null) sources.add(alphabet);
      if (posSeq != null) sources.add(posSeq);
      if (negSeq != null) sources.add(negSeq);
      return sources;
    }
  
    @Override
    public void cleanUp() {
      if (alphabet != null) alphabet.getFile().delete();
      if (posSeq != null) posSeq.getFile().delete();
      if (negSeq != null) negSeq.getFile().delete();
    }
  }

  public Dreme() {
    super("DREME", "DREME");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the template
    this.tmplMain = cache.loadAndCache("/WEB-INF/templates/dreme.tmpl");
    this.tmplVerify = cache.loadAndCache("/WEB-INF/templates/dreme_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    control = new ComponentSequences(context, tmplMain.getSubtemplate("control"));
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
        jobDetails.getHelp(), advancedOptions.getHelp(), submitReset.getHelp(),
        footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("alphabet", alphabet.getComponent());
    main.set("sequences", sequences.getComponent());
    main.set("control", control.getComponent());
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
      // stopping criterion
      if ("enable_nmotifs".equals(paramRequire(request, "srch_limit"))) {
        // motif count
        data.nMotifs = paramInteger(feedback, "Motif count", request, "nmotifs", 1, null, 10);
        data.ethresh = 0;
      } else {
        // E-value threshold
        data.ethresh = paramNumber(feedback, "<i>E</i>-value threshold", request,
          "ethresh", 0.0, null, 0.05);
      }
      // no reverse complement
      data.norc = paramBool(request, "norc");
      error = false;
    } finally {
      if (error) data.cleanUp();
    }
    return data;
  }

}

