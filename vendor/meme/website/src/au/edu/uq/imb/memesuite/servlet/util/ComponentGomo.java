package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import java.sql.SQLException;
import java.util.EnumSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.paramRequire;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.GOMO_DB_KEY;

/**
 * Data entry component for GOMo databases
 */
public class ComponentGomo extends PageComponent {
  private ServletContext context;
  private HTMLTemplate tmplGomo;
  private HTMLTemplate tmplCategory;
  private HTMLTemplate tmplListing;
  private String prefix;
  private HTMLTemplate title;
  private HTMLTemplate subtitle;


  private static final Pattern DB_ID_PATTERN = Pattern.compile("^\\d+$");
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.gomo");

  public ComponentGomo(ServletContext context, HTMLTemplate info) throws ServletException {
    this.context = context;
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplGomo = cache.loadAndCache("/WEB-INF/templates/component_gomo.tmpl");
    prefix = getText(info, "prefix", "motifs");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    tmplCategory = tmplGomo.getSubtemplate("component").getSubtemplate("category");
    tmplListing = tmplCategory.getSubtemplate("listing");
  }

  @Override
  public HTMLSub getComponent() {
    HTMLSub gomo = tmplGomo.getSubtemplate("component").toSub();
    gomo.set("prefix", prefix);
    if (title != null) gomo.set("title", title);
    if (subtitle != null) gomo.set("subtitle", subtitle);
    GomoDBList db = (GomoDBList)context.getAttribute(GOMO_DB_KEY);
    if (db == null) {
      logger.log(Level.SEVERE, "No GOMo database is intalled");
      gomo.empty("category");
      return gomo;
    }
    try {
      gomo.set("category", new AllCategories(db));
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Error querying GOMo categories", e);
      gomo.empty("category");
    }
    return gomo;
  }

  @Override
  public HTMLSub getHelp() {
    return this.tmplGomo.getSubtemplate("help").toSub();
  }

  public GomoDB getGomoSequences(HttpServletRequest request) throws ServletException {
    // determine the source
    String source = paramRequire(request, prefix + "_source");
    Matcher m = DB_ID_PATTERN.matcher(source);
    if (!m.matches()) {
      throw new ServletException("Parameter " + prefix + "_source had a " +
          "value that did not match any of the allowed values.");
    }
    long dbId;
    try {
      dbId = Long.parseLong(source, 10);
    } catch (NumberFormatException e) {
      throw new ServletException(e);
    }
    GomoDBList db = (GomoDBList)context.getAttribute(GOMO_DB_KEY);
    if (db == null) {
      throw new ServletException("Unable to access the gomo database.");
    }
    try {
      return db.getGomoListing(dbId);
    } catch (SQLException e) {
      throw new ServletException(e);
    }
  }

  private class AllCategories extends HTMLSubGenerator<Category> {
    private GomoDBList db;
    private AllCategories(GomoDBList db) throws SQLException {
      super(db.getCategories(false, EnumSet.allOf(AlphStd.class)));
      this.db = db;
    }
    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = tmplCategory.toSub();
      out.set("name", item.getName());
      try {
        out.set("listing", new AllListingsOfCategory(db, item.getID()));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error querying GOMo listings", e);
        out.empty("listing");
      }
      return out;
    }
  }

  private class AllListingsOfCategory extends HTMLSubGenerator<Listing> {
    private AllListingsOfCategory(GomoDBList db, long id) throws SQLException {
      super(db.getListings(id));
    }

    @Override
    protected HTMLSub transform(Listing item) {
      HTMLSub out = tmplListing.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      return out;
    }
  }
}
