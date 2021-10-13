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
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.TGENE_DB_KEY;

/**
 * Data entry component for Tgene databases
 */
public class ComponentTgene extends PageComponent {
  private ServletContext context;
  private HTMLTemplate tmplTgene;
  private HTMLTemplate tmplCategory;
  private HTMLTemplate tmplListing;
  private String prefix;
  private HTMLTemplate title;
  private HTMLTemplate subtitle;


  private static final Pattern DB_ID_PATTERN = Pattern.compile("^\\d+$");
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.tgene");

  public ComponentTgene(ServletContext context, HTMLTemplate info) throws ServletException {
    this.context = context;
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplTgene = cache.loadAndCache("/WEB-INF/templates/component_tgene.tmpl");
    prefix = getText(info, "prefix", "motifs");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    tmplCategory = tmplTgene.getSubtemplate("component").getSubtemplate("category");
    tmplListing = tmplCategory.getSubtemplate("listing");
  }

  @Override
  public HTMLSub getComponent() {
    HTMLSub tgene = tmplTgene.getSubtemplate("component").toSub();
    tgene.set("prefix", prefix);
    if (title != null) tgene.set("title", title);
    if (subtitle != null) tgene.set("subtitle", subtitle);
    TgeneDBList db = (TgeneDBList)context.getAttribute(TGENE_DB_KEY);
    if (db == null) {
      logger.log(Level.SEVERE, "No Tgene database is intalled");
      tgene.empty("category");
      return tgene;
    }
    try {
      tgene.set("category", new AllCategories(db));
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Error querying Tgene categories", e);
      tgene.empty("category");
    }
    return tgene;
  }

  @Override
  public HTMLSub getHelp() {
    return this.tmplTgene.getSubtemplate("help").toSub();
  }

  public TgeneDB getTgenePanel(HttpServletRequest request) throws ServletException {
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
    TgeneDBList db = (TgeneDBList)context.getAttribute(TGENE_DB_KEY);
    if (db == null) {
      throw new ServletException("Unable to access the Tgene database.");
    }
    try {
      return db.getTgeneListing(dbId);
    } catch (SQLException e) {
      throw new ServletException(e);
    }
  } // getTgenePanel

  private class AllCategories extends HTMLSubGenerator<Category> {
    private TgeneDBList db;
    private AllCategories(TgeneDBList db) throws SQLException {
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
        logger.log(Level.SEVERE, "Error querying Tgene listings", e);
        out.empty("listing");
      }
      return out;
    }
  }

  private class AllListingsOfCategory extends HTMLSubGenerator<Listing> {
    private AllListingsOfCategory(TgeneDBList db, long id) throws SQLException {
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
} //class ComponentTgene
