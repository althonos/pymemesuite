package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.servlet.util.WebUtils;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

// debugging
//import java.lang.Object;
//import java.io.Serializable;
//import org.apache.commons.lang.builder.ToStringBuilder;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.TGENE_DB_KEY;

/**
 * Show the databases available to the T-Gene program.
 */
public class ShowTgeneDBs extends HttpServlet {
  private ServletContext context;
  private HTMLTemplate template;
  private HTMLTemplate categoryTemplate;
  private HTMLTemplate listingTemplate;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web");

  public ShowTgeneDBs() { }

  @Override
  public void init() throws ServletException {
    context = this.getServletContext();
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/show_tgene_dbs.tmpl");
    categoryTemplate = template.getSubtemplate("content").getSubtemplate("category");
    listingTemplate = categoryTemplate.getSubtemplate("listing");
  }

  @Override
  public void doGet(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    if (request.getParameter("category") != null) {
      outputXmlListingsOfCategory(response, getId(request, "category"));
    } else {
      display(response);
    }
  }

  @Override
  public void doPost(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    display(response);
  }

  private void display(HttpServletResponse response)
      throws IOException, ServletException {
    TgeneDBList tgeneDBList = (TgeneDBList)context.getAttribute(TGENE_DB_KEY);
    response.setContentType("text/html; charset=UTF-8");
    HTMLSub out = template.toSub();
    if (tgeneDBList != null) {
      try {
        out.getSub("content").set("category", new AllCategories(tgeneDBList));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading categories", e);
        out.empty("content");
      }
    } else {
      out.empty("content");
    }
    out.output(response.getWriter());
  }

  private long getId(HttpServletRequest request, String name) throws ServletException {
    String value = request.getParameter(name);
    if (value == null) {
      throw new ServletException("Parameter '" + name + "' was not set.");
    }
    long id;
    try {
      id = Long.parseLong(value, 10);
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter '" + name + "' is not a integer value.", e);
    }
    return id;
  }

  private void outputXmlListingsOfCategory(HttpServletResponse response, long categoryId)
      throws IOException, ServletException {
    TgeneDBList tgeneDBList = (TgeneDBList)context.getAttribute(TGENE_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    PrintWriter out = null;
    try {
      Iterator<Listing> listingView = tgeneDBList.getListings(categoryId).iterator();
      out = response.getWriter();
      out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      out.println("<listings category=\"" + categoryId + "\">");
      while (listingView.hasNext()) {
        Listing listing = listingView.next();
        out.println("<l i=\"" + listing.getID() + "\" n=\"" + WebUtils.escapeForXML(listing.getName()) + "\"/>");
      }
      out.println("</listings>");
    } catch (SQLException e) {
      throw new ServletException("Error getting T-Gene listings", e);
    } finally {
      if (out != null) out.close();
    }
  }

  private class AllCategories extends HTMLSubGenerator<Category> {
    private TgeneDBList tgeneDBList;
    private AllCategories(TgeneDBList tgeneDBList) throws SQLException {
      super(tgeneDBList.getCategories(false, EnumSet.allOf(AlphStd.class)));
      this.tgeneDBList = tgeneDBList;
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      try {
        out.set("listing", new AllListingsOfCategory(tgeneDBList, item.getID()));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading T-Gene listings", e);
        out.empty("listing");
      }
      return out;
    }
  }

  private class AllListingsOfCategory extends HTMLSubGenerator<Listing> {
    private TgeneDBList tgeneDBList;
    private AllListingsOfCategory(TgeneDBList tgeneDBList, long id) throws SQLException{
      super(tgeneDBList.getListings(id));
      this.tgeneDBList = tgeneDBList;
    }

    @Override
    protected HTMLSub transform(Listing listing) {
      HTMLSub out = listingTemplate.toSub();
      out.set("name", listing.getName());
      out.set("description", listing.getDescription());
      try {
        TgeneDB tgeneDB = tgeneDBList.getTgeneListing(listing.getID());
        boolean genome_only = tgeneDB.getTissues().equals("none");
        out.set("genome_release", tgeneDB.getGenomeRelease());
        out.set("annotation_file", tgeneDB.getAnnotationFileName());
        out.set("tissues", genome_only ? null : tgeneDB.getTissues().replace(",", " "));
        out.set("histones", tgeneDB.getHistones().replace(",", " "));
        out.set("max_link_distances", tgeneDB.getMaxLinkDistances().replace(",", " "));
        out.set("rna_source", tgeneDB.getRnaSource());
        out.set("lecat", tgeneDB.getLecat());
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error querying T-Gene listing",  e);
      }
      return out;
    }
  }

}
