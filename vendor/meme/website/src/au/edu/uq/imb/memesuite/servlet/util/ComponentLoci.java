package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.*;
//import au.edu.uq.imb.memesuite.io.bed.BedException;
//import au.edu.uq.imb.memesuite.io.bed.BedParser;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.Part;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;


/**
 * A component that is used for inputting chromosome loci in BED format.
 */
public class ComponentLoci extends PageComponent {
  private ServletContext context;
  private HTMLTemplate tmplLoci;
  private String prefix;
  private String fieldName; // for feedback
  private HTMLTemplate title;
  private HTMLTemplate subtitle;
  private boolean enableNoLoci;
  private String registerFn;
  private DefaultOption defaultOption;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.locus");

  /**
   * Enum of the initial selection state.
   */
  public static enum DefaultOption {
    NOLOCI,
    TEXT,
    FILE,
  }

  public ComponentLoci(ServletContext context, HTMLTemplate info) throws ServletException {
    this.context = context;
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplLoci = cache.loadAndCache("/WEB-INF/templates/component_loci.tmpl");
    prefix = getText(info, "prefix", "loci");
    fieldName = getText(info, "description", "loci");
    registerFn = getText(info, "register", "nop");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    enableNoLoci = info.containsSubtemplate("enable_noloci");
  } // ComponentLoci

  public HTMLSub getComponent() {
    DefaultOption defaultOption = this.defaultOption;
    HTMLSub loci = tmplLoci.getSubtemplate("component").toSub();
    loci.set("prefix", prefix);
    if (title != null) loci.set("title", title);
    if (subtitle != null) loci.set("subtitle", subtitle);

    StringWriter buf = new StringWriter();
    JsonWr jsonWr = new JsonWr(buf, 18);
    try {
      jsonWr.start();
      jsonWr.property("field", fieldName);
      jsonWr.end();
    } catch (IOException e) {
      // no IO exceptions should occur as this uses a StringBuffer
      throw new Error(e);
    }
    loci.set("register_component", registerFn);
    return loci;
  }

  public HTMLSub getHelp() {
    return tmplLoci.getSubtemplate("help").toSub();
  }

  private boolean checkSpec(LociStats stats, FeedbackHandler feedback) {
    boolean ok = true;
    return ok;
  }

  /**
   * Loci come from 1 sources.
   * They can be uploaded.
   * Uploaded loci need to be preprocessed to calculate
   * statistics and to ensure they are valid.
   *
   * So what should happen if...
   * 1) An expected field is missing
   *    throw an exception and stop processing the request
   * 2) Loci have not been sent
   *    whine complaining about the missing loci
   *    return null
   * 3) Loci contain a syntax error
   *    whine complaining about the syntax error
   *    return null
   * 4) Loci violate some constraint
   *    whine complaining about the failed constraint
   *    return null
   * 5) Loci are valid, pass all constraints
   *    return a locus source
   *
   *
   * @param name manages the name of user supplied loci and avoids clashes.
   * @param request all the information sent to the webserver.
   * @param feedback an interface for providing error messages to the user.
   * @return a locus source and return null when a source is not available.
   * @throws ServletException if request details are incorrect (like missing form fields).
   * @throws IOException if storing a parsed version of the loci to file fails.
   */
  public LociInfo getLoci(FileCoord.Name name,
      HttpServletRequest request, FeedbackHandler feedback) throws ServletException, IOException {
      Part part = request.getPart(prefix + "_file");
      if (part == null || part.getSize() == 0) {
        feedback.whine("No " + fieldName + " provided.");
        return null; // no loci submitted
      }
      name.setOriginalName(getPartFilename(part));
      LociStats statistics = null;
      InputStream in = null;
      File file = null;
      OutputStream out = null;
      boolean success = false;
      byte[] buffer = new byte[1000];
      try {
        in  = new BufferedInputStream(part.getInputStream());
        file = File.createTempFile("uploaded_loci_", ".bed");
        file.deleteOnExit();
        out = new BufferedOutputStream(new FileOutputStream(file));
        //BedParser.parse(in, handler);
	// Just copy the input BED file for now; FIXME
        int numBytes;
	while ((numBytes = in.read(buffer)) != -1) {
          out.write(buffer, 0, numBytes);
        }
        try {in.close();} finally {in = null;}
        try {out.close();} finally {out = null;}
        // check the statistics
        //statistics = handler.getStatsRecorder();
        statistics = new LociStats();
        //success = checkSpec(statistics, feedback);
        success = true;
      //} catch (BedException e) {
        // ignore
      } finally {
        closeQuietly(in);
        if (file != null && !success) {
          if (!file.delete()) file.deleteOnExit();
        }
      }
      if (success) return new LociDataSource(file, name, statistics);
      return null;
    //}

  } // getLoci

} // class ComponentLoci
