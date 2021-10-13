package au.edu.uq.imb.memesuite.db;

import java.io.*;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.*;
import java.sql.*;

import au.edu.uq.imb.memesuite.data.AlphStd;
import org.sqlite.*;

public class TgeneDBList extends DBList {
  private File tsv;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.tgenedb");

  public TgeneDBList(File tsv) throws ClassNotFoundException, IOException, SQLException {
    super(loadTgeneTSV(tsv), true);
    this.tsv = tsv;
  }

  /**
   * Get the source file that the database was created from.
   * @return the TSV file that was loaded to create the database.
   */
  public File getTSV() {
    return tsv;
  }

  @Override
  protected PreparedStatement prepareListingsQuery(Connection conn, long categoryId, boolean shortOnly, EnumSet<AlphStd> allowedAlphabets) throws SQLException {
    PreparedStatement ps = conn.prepareStatement(SQL.SELECT_TGENE_LISTINGS_OF_CATEGORY);
    ps.setLong(1, categoryId);
    return ps;
  }

  public TgeneDB getTgeneListing(long listingId) throws SQLException {
    TgeneDB tgeneDB;
    Connection conn = null;
    PreparedStatement stmt = null;
    ResultSet rset = null;
    try {
      // create a database connection
      int i = 1;
      conn = ds.getConnection();
      stmt = conn.prepareStatement(SQL.SELECT_TGENE_LISTING);
      stmt.setLong(1, listingId);
      rset = stmt.executeQuery();
      if (!rset.next()) throw new SQLException("No listing by that ID");
      String name = rset.getString(i++);
      String description = rset.getString(i++);
      String genome_release = rset.getString(i++);
      String rna_source = rset.getString(i++);
      String tissues = rset.getString(i++);
      String histone_root = rset.getString(i++);
      String histones = rset.getString(i++);
      String max_link_distances = rset.getString(i++);
      String expression_root = rset.getString(i++);
      String annotation_file_name = rset.getString(i++);
      String transcript_types = rset.getString(i++);
      String use_gene_ids = rset.getString(i++);
      double lecat = rset.getDouble(i++);
      rset.close(); rset = null;
      stmt.close(); stmt = null;
      conn.close(); conn = null;
      tgeneDB = new TgeneDB( listingId, name, description, genome_release, 
	rna_source, tissues, histone_root, histones, max_link_distances, 
	expression_root, 
        annotation_file_name, 
        transcript_types, use_gene_ids, lecat
      );
    } finally {
      if (rset != null) {
        try {
          rset.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (stmt != null) {
        try {
          stmt.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (conn != null) {
        try {
          conn.close();
        } catch (SQLException e) {
          logger.log(Level.SEVERE, "Cleanup (after error) failed to close the DB connection", e);
        }
      }
    }
    return tgeneDB;
  } // getTgeneListing

  private static File loadTgeneTSV(File tsv) throws ClassNotFoundException, IOException, SQLException {
    // load the JDBC needed
    Class.forName("org.sqlite.JDBC");
    // create the file to contain the database
    File db = File.createTempFile("tgene_db_", ".sqlite");
    db.deleteOnExit();
    // configure the database
    SQLiteConfig config = new SQLiteConfig();
    config.enforceForeignKeys(true);
    SQLiteDataSource ds = new SQLiteDataSource(config);
    ds.setUrl("jdbc:sqlite:" + db);
    // open a connection
    Connection connection = null;
    try {
      connection = ds.getConnection();
      Statement statement = connection.createStatement();
      statement.executeUpdate(SQL.CREATE_TBL_CATEGORY);
      statement.executeUpdate(SQL.CREATE_TBL_LISTING);
      statement.executeUpdate(SQL.CREATE_TBL_TGENE);
      connection.setAutoCommit(false);
      importTgeneTSV(tsv, connection);
      connection.commit();
    } finally {
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    if (!db.setLastModified(tsv.lastModified())) {
      logger.log(Level.WARNING, "Failed to set last modified date on " + db);
    }
    logger.log(Level.INFO, "Loaded TGENE TSV \"" + tsv + "\" into \"" + db + "\"");
    return db;
  } // loadTgeneTSV

  private static void importTgeneTSV(File tsv, Connection connection) throws IOException, SQLException {
    String line;
    Long categoryId = null;
    Pattern emptyLine = Pattern.compile("^\\s*$");
    Pattern commentLine = Pattern.compile("^#");
    Pattern dashedName = Pattern.compile("^-*([^-](?:.*[^-])?)-*$");
    BufferedReader in = null;
    try {
      ResultSet rs;
      // open the tsv file for reading
      in = new BufferedReader(new InputStreamReader(new FileInputStream(tsv), "UTF-8"));
      // create the prepared statements
      PreparedStatement pstmtCategory = connection.prepareStatement(
          SQL.INSERT_CATEGORY, Statement.RETURN_GENERATED_KEYS);
      PreparedStatement pstmtListing = connection.prepareStatement(
          SQL.INSERT_LISTING, Statement.RETURN_GENERATED_KEYS);
      PreparedStatement pstmtTgeneDB = connection.prepareStatement(
          SQL.INSERT_TGENE_DB, Statement.NO_GENERATED_KEYS);
      // now read the tsv file
      while ((line = in.readLine()) != null) {
        // skip any empty lines or comments
        if (emptyLine.matcher(line).find()) continue;
        if (commentLine.matcher(line).find()) continue;
        String[] values = line.split("\\t");	// split on tab
        // check that a name was supplied
        if (emptyLine.matcher(values[0]).find()) {
          logger.log(Level.WARNING, "TSV line has no entry for name column.");
          continue;
        }
        // test to see if we have a category or a selectable listing
        if (values.length == 1) {
          // category
          String name = values[0].trim();
          // remove dashes from around name
          Matcher matcher = dashedName.matcher(name);
          if (matcher.matches()) {
            name = matcher.group(1);
          }
          pstmtCategory.setString(1, name);
          pstmtCategory.executeUpdate();
          // now get the category Id
          rs = pstmtCategory.getGeneratedKeys();
          if (!rs.next()) throw new SQLException("Failed to get Category Id.\n");
          categoryId = rs.getLong(1);
          rs.close();
          logger.log(Level.FINE, "Loaded T-Gene Panel Category: " + name);
        } else {
	  // listing
	  if (values.length < 13) {
	    logger.log(Level.WARNING, "TSV line has " + values.length +
		" values but expected 13 for a TGENE listing.");
            logger.log(Level.WARNING, "Line is: '" + line + "'");
	    continue;
          }
	  if (categoryId == null) {
	    // create a dummy category with an empty name
	    pstmtCategory.setString(1, "");
	    pstmtCategory.executeUpdate();
	    // now get the category Id
	    rs = pstmtCategory.getGeneratedKeys();
	    if (!rs.next()) throw new SQLException("Failed to get Category Id.\n");
	    categoryId = rs.getLong(1);
	    rs.close();
	  }
	  int i = 0;
	  String name = values[i++];
	  String description = values[i++];
	  String genome_release = values[i++];
	  String rna_source = values[i++];
	  String tissues = values[i++];
	  String histone_root = values[i++];
	  String histones = values[i++];
	  String max_link_distances = values[i++];
	  String expression_root = values[i++];
	  String annotation_file_name = values[i++];
	  String transcript_types = values[i++];
	  String use_gene_ids = values[i++];
	  double lecat = Double.parseDouble(values[i++]);
	  // check we have a category to store the listing under
	  // now create the listing
	  pstmtListing.setLong(1, categoryId);
	  pstmtListing.setString(2, name);
	  pstmtListing.setString(3, description);
	  pstmtListing.executeUpdate();
	  rs = pstmtListing.getGeneratedKeys();
	  if (!rs.next()) throw new SQLException("Failed to get Listing Id.\n");
	  long listingId = rs.getLong(1);
	  rs.close();
	  // create the TGENE primary listing
	  i = 1; 
	  pstmtTgeneDB.setLong(i++, listingId);
	  pstmtTgeneDB.setString(i++, genome_release);
	  pstmtTgeneDB.setString(i++, rna_source);
	  pstmtTgeneDB.setString(i++, tissues);
	  pstmtTgeneDB.setString(i++, histone_root);
	  pstmtTgeneDB.setString(i++, histones);
	  pstmtTgeneDB.setString(i++, max_link_distances);
	  pstmtTgeneDB.setString(i++, expression_root);
	  pstmtTgeneDB.setString(i++, annotation_file_name);
	  pstmtTgeneDB.setString(i++, transcript_types);
	  pstmtTgeneDB.setString(i++, use_gene_ids);
	  pstmtTgeneDB.setDouble(i++, lecat);
	  pstmtTgeneDB.executeUpdate();
	  logger.log(Level.FINE, "Loaded Tgene Listing: " + name);
        }
      }
      // close the prepared statements
      pstmtCategory.close();
      pstmtListing.close();
      pstmtTgeneDB.close();
    } finally {
      if (in != null) {
        try {
          in.close();
        } catch (IOException e) {
          // ignore
        }
      }
    } // try
  } // importTgeneTSV
} // class TgeneDBList 
