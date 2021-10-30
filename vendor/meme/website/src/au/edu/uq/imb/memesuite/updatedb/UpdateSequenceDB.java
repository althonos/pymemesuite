package au.edu.uq.imb.memesuite.updatedb;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.util.GlobFilter;
import au.edu.uq.imb.memesuite.util.MultiSourceStatus;
import au.edu.uq.imb.memesuite.util.Progress;
import com.martiansoftware.jsap.*;
import com.martiansoftware.jsap.stringparsers.FileStringParser;
import org.sqlite.SQLiteConfig;
import org.sqlite.SQLiteDataSource;

import java.io.*;
import java.lang.reflect.Constructor;
import java.sql.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Date;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static au.edu.uq.imb.memesuite.db.SQL.*;
import static au.edu.uq.imb.memesuite.util.MultiSourceStatus.*;

/**
 * Updates sequence databases
 */
public class UpdateSequenceDB {

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.updatedb");

  /**
   * Creates a SQLite DataSource for the file.
   * @param dbFile the SQLite database file.
   * @return the created DataSource.
   * @throws ClassNotFoundException when the SQLite library cannot be found in the classpath
   */
  private static SQLiteDataSource getDBSource(File dbFile) throws ClassNotFoundException {
    SQLiteDataSource ds;
    Class.forName("org.sqlite.JDBC");
    SQLiteConfig config = new SQLiteConfig();
    config.enforceForeignKeys(true);
    config.setReadOnly(false);
    ds = new SQLiteDataSource(config);
    ds.setUrl("jdbc:sqlite:" + dbFile);
    return ds;
  }

  /**
   * Creates the tables required to store information on the downloaded sequence
   * databases. This assumes a freshly created database.
   * @param ds the DataSource.
   * @throws SQLException when the SQL doesn't work.
   */
  private static void setupTables(SQLiteDataSource ds) throws SQLException {
    logger.log(Level.INFO, "Creating database tables.");
    Connection connection = null;
    try {
      connection = ds.getConnection();
      Statement statement = connection.createStatement();
      logger.log(Level.INFO, "Creating category table");
      statement.executeUpdate(CREATE_TBL_CATEGORY);
      logger.log(Level.INFO, "Creating listing table");
      statement.executeUpdate(CREATE_TBL_LISTING);
      logger.log(Level.INFO, "Creating sequence table");
      statement.executeUpdate(CREATE_TBL_SEQUENCE_FILE);
      logger.log(Level.INFO, "Creating priors table");
      statement.executeUpdate(CREATE_TBL_PRIOR_FILE);
      logger.log(Level.INFO, "Declaring the schema version");
      statement.executeUpdate(DECLARE_VERSION);
      statement.close();
      connection.close();
      connection = null;
    } finally {
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) {
          logger.log(Level.WARNING, "Failed to close connection to database.");
        }
      }
    }
    logger.log(Level.INFO, "Successfully created database tables.");
  }

  /**
   * Allows testing a table schema to validate a database table.
   */
  private static class TableSchema {
    private String name;
    private Map<String,String> columns;

    /**
     * Construct a schema for the named table.
     * @param name the name of the table.
     */
    public TableSchema(String name) {
      this.name = name;
      columns = new TreeMap<String, String>();
    }

    /**
     * Add a checked column to this table schema
     * @param name the column name.
     * @param type the column type.
     * @return this TableSchema so the method can be chained.
     */
    public TableSchema addCol(String name, String type) {
      columns.put(name, type);
      return this;
    }

    /**
     * Checks this table schema against the SQLite database this is connected to.
     * @param connection a connection to a SQLite database (uses SQLite specific query)
     * @return true if the table schema matches the table in the connected SQLite database.
     * @throws SQLException if the query required to check the schema fails.
     */
    public boolean validate(Connection connection) throws SQLException {
      Statement statement = connection.createStatement();
      ResultSet resultSet = statement.executeQuery("PRAGMA table_info(" + this.name + ")");
      int columnCount = 0;
      boolean ok = true;
      while (resultSet.next()) {
        String col_name = resultSet.getString(2);
        String col_type = resultSet.getString(3);
        if (columns.containsKey(col_name)) {
          if (!columns.get(col_name).equals(col_type)) {
            ok = false;
            logger.log(Level.WARNING, "Table " + name + " column " + col_name +
                " has unexpected type. Expected " + columns.get(col_name) + " but got " + col_type);
          }
          columnCount++;
        }
      }
      if (columnCount < this.columns.size()) {
        logger.log(Level.WARNING, "Table " + name + " is missing expected columns.");
        ok = false;
      }
      return ok;
    }
  }

  /**
   * Reads the version of the database schema.
   * @param connection a connection to a SQLite database (uses SQLite specific query).
   * @return the version or 0 if it is missing.
   * @throws SQLException if an error occurs.
   */
  private static int queryDbVersion(Connection connection) throws SQLException {
    int version = 0;
    Statement statement = connection.createStatement();
    ResultSet resultSet = statement.executeQuery(QUERY_VERSION);
    if (resultSet.next()) {
      version = resultSet.getInt(1);
    }
    resultSet.close();
    statement.close();
    return version;
  }

  private static void upgradeV1ToV2(SQLiteDataSource ds) throws SQLException {
    logger.log(Level.INFO, "Updating table schema from version 1 to version 2");
     // need to add the obsolete column
    Connection conn = null;
    try {
      conn = ds.getConnection();
      conn.setAutoCommit(false);
      // alter tblSequenceFile to add obsolete column
      Statement stmt = conn.createStatement();
      stmt.executeUpdate("ALTER TABLE tblSequenceFile ADD COLUMN obsolete INTEGER DEFAULT 0");
      stmt.executeUpdate(SQL.DECLARE_VERSION);
      stmt.close();
      conn.commit();
      conn.close();
      conn = null;
      logger.log(Level.INFO, "Successfully updated table schema from version 1 to version 2");
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Failed to update table schema from version 1 to version 2");
      throw e;
    } finally {
      if (conn != null) {
        try {
          conn.rollback();
        } catch (SQLException e) {/* ignore */}
        try {
          conn.close();
        } catch (SQLException e) {/* ignore */}
      }
    }
  }

  private static void upgradeV2ToV3(SQLiteDataSource ds) throws SQLException {
    logger.log(Level.INFO, "Updating table schema from version 2 to version 3");
    // need to add the obsolete column
    Connection conn = null;
    try {
      conn = ds.getConnection();
      conn.setAutoCommit(false);
      // alter tblSequenceFile to add obsolete column
      Statement stmt = conn.createStatement();
      stmt.executeUpdate(SQL.CREATE_TBL_PRIOR_FILE);
      stmt.executeUpdate(SQL.DECLARE_VERSION);
      stmt.close();
      conn.commit();
      conn.close();
      conn = null;
      logger.log(Level.INFO, "Successfully updated table schema from version 2 to version 3");
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Failed to update table schema from version 2 to version 3");
      throw e;
    } finally {
      if (conn != null) {
        try {
          conn.rollback();
        } catch (SQLException e) {/* ignore */}
        try {
          conn.close();
        } catch (SQLException e) {/* ignore */}
      }
    }
  }

  /**
   * Checks that the database schema is valid for a sequence database.
   * @param ds the database data source.
   * @return true if the version number is compatible and all the expected tables are present with all the expected columns.
   * @throws SQLException if a query goes wrong.
   */
  private static boolean validateTables(SQLiteDataSource ds) throws SQLException {
    boolean ok = true;
    TableSchema[] schemas = new TableSchema[] {
        new TableSchema("tblCategory").addCol("id", "INTEGER").addCol("name", "TEXT"),
        new TableSchema("tblListing").addCol("id", "INTEGER").
            addCol("categoryId", "INTEGER").addCol("name", "TEXT").
            addCol("description", "TEXT"),
        new TableSchema("tblSequenceFile").addCol("id", "INTEGER").
            addCol("retriever", "INTEGER").addCol("listingId", "INTEGER").
            addCol("alphabet", "INTEGER").addCol("edition", "INTEGER").
            addCol("version", "TEXT").addCol("description", "TEXT").
            addCol("fileSeq", "TEXT").addCol("fileBg", "TEXT").
            addCol("sequenceCount", "INTEGER").addCol("totalLen", "INTEGER").
            addCol("minLen", "INTEGER").addCol("maxLen", "INTEGER").
            addCol("avgLen", "REAL").addCol("stdDLen", "REAL").
            addCol("obsolete", "INTEGER"),
        new TableSchema("tblPriorFile").addCol("id", "INTEGER").
            addCol("sequenceId", "INTEGER").addCol("biosample", "TEXT").
            addCol("assay", "TEXT").addCol("source", "TEXT").
            addCol("description", "TEXT").addCol("filePrior", "TEXT")
    };
    Connection connection = null;
    try {
      connection = ds.getConnection();
      int version = queryDbVersion(connection);
      if (version == 0) {
        logger.log(Level.SEVERE, "The database schema version is missing.");
        ok = false;
      } else if (version < SCHEMA_VERSION) {
        switch (version) {
          case 1:
            upgradeV1ToV2(ds);
          case 2:
            upgradeV2ToV3(ds);
          case 3:
            break;
          default:
            logger.log(Level.SEVERE, "The database schema is older than handled by this program.");
            ok = false;
        }
      } else if (version > SCHEMA_VERSION) {
        logger.log(Level.SEVERE, "The database schema is newer than handled by this program.");
        ok = false;
      }
      // check that at least the fields we require are present
      for (TableSchema schema : schemas) {
        ok &= schema.validate(connection);
      }
      connection.close();
      connection = null;
    } finally {
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) {
          logger.log(Level.WARNING, "Failed to close connection to database.");
        }
      }
    }
    return ok;
  }

  /**
   * Gets the database source when given the database directory.
   * @param dbDir the directory containing the database
   * @return a database data source.
   * @throws SQLException when queries fail.
   * @throws IOException when files can't be created/read.
   * @throws ClassNotFoundException when the database library can't be loaded.
   */
  private static SQLiteDataSource getInitialisedDatabase(File dbDir)
      throws SQLException, IOException, ClassNotFoundException {
    boolean newDB;
    File db = new File(dbDir, "fasta_db.sqlite");
    if ((newDB = db.createNewFile())) {
      logger.log(Level.INFO, "Creating a new sqlite database " + db);
    }
    SQLiteDataSource ds = getDBSource(db);
    if (newDB) {
      setupTables(ds);
    } else {
      if (!validateTables(ds)) {
        throw new SQLException("Incompatible or invalid database.");
      }
    }
    return ds;
  }

  public static void loadPriors(SQLiteDataSource ds, File dbDir, File priorsTsv, MultiSourceStatus status) {
    Connection conn = null;
    BufferedReader reader = null;
    Pattern splitter = Pattern.compile("\t");
    Pattern entry = Pattern.compile("\\S");
    Progress progress = new Progress(status, "Loading Priors");
    try {
      conn = ds.getConnection();
      conn.setAutoCommit(false);
      // Drop and recreate priors table
      Statement stmt = conn.createStatement();
      stmt.executeUpdate(SQL.DROP_TBL_PRIOR_FILE);
      stmt.executeUpdate(SQL.CREATE_TBL_PRIOR_FILE);
      stmt.close();
      // Prepare to load the prior files
      PreparedStatement querySeqId = conn.prepareStatement(SQL.SELECT_SEQUENCE_FILE_ID);
      PreparedStatement insertPriorFile = conn.prepareStatement(SQL.INSERT_PRIOR_FILE);
      // open the file
      reader = new BufferedReader(new InputStreamReader(new FileInputStream(priorsTsv)));
      String line;
      int lineNo = 0;
      while ((line = reader.readLine()) != null) {
        logger.log(Level.INFO, "Read line\n" + line);
        lineNo++;
        // skip empty lines
        if (!entry.matcher(line).find()) {
          logger.log(Level.INFO, "Skipping empty line");
          continue;
        }
        // skip comments
        if (line.charAt(0) == '#') {
          logger.log(Level.INFO, "Skipping comment line");
          continue;
        }
        // split into fields
        String[] fields = splitter.split(line, -1);
        // skip when there are not enough fields or too many fields
        if (fields.length != 8) {
          logger.log(Level.WARNING, "Skipping line " + lineNo + " in priors configuration because it does not list exactly 8 fields");
          continue;
        }
        // assign the fields to nicer names
        String sequenceFileName = fields[0];
        String priorsFileName = fields[1];
        String priorsDistFileName = fields[2];
        String biosample = fields[3];
        String assay = fields[4];
        String source = fields[5];
        String url = fields[6];
        String description = fields[7];
        // check for existence of the priors file
        File priorsFile = new File(dbDir, priorsFileName);
        if (!priorsFile.exists()) {
          logger.log(Level.WARNING, "Skipping line " + lineNo + " in priors configuration because the priors file \"" +
              priorsFileName + "\" does not exist in the database directory");
          continue;
        }
        // check for existence of the priors-dist file
        File priorsDistFile = new File(dbDir, priorsDistFileName);
        if (!priorsDistFile.exists()) {
          logger.log(Level.WARNING, "Skipping line " + lineNo + " in priors configuration because the priors-dist file \"" +
              priorsDistFileName + "\" does not exist in the database directory");
          continue;
        }
        // lookup the sequence ID
        querySeqId.setString(1, sequenceFileName);
        ResultSet rset = querySeqId.executeQuery();
        if (!rset.next()) {
          logger.log(Level.INFO, "Skipping line " + lineNo + " in priors configuration because the sequence file \"" +
              sequenceFileName + "\" is not listed in the database");
          rset.close();
          continue;
        }
        long sequenceId = rset.getLong(1);
        rset.close();
        // insert a new entry into the database
        insertPriorFile.setLong(1, sequenceId);
        insertPriorFile.setString(2, priorsFileName);
        insertPriorFile.setString(3, priorsDistFileName);
        insertPriorFile.setString(4, biosample);
        insertPriorFile.setString(5, assay);
        insertPriorFile.setString(6, source);
        insertPriorFile.setString(7, url);
        insertPriorFile.setString(8, description);
        insertPriorFile.execute();
        logger.log(Level.INFO, "Created prior for sequence file " + sequenceFileName);
      }
      reader.close();
      reader = null;
      // close statements
      insertPriorFile.close();
      querySeqId.close();
      // all done, commit changes
      progress.setTask("Commiting changes");
      conn.commit();
      conn.close();
      conn = null;
      logger.log(Level.INFO, "Finished loading priors.");
    } catch (SQLException e) {
      progress.setTask("Error - check logs!");
      logger.log(Level.SEVERE, "Failed to load priors", e);
    } catch (IOException e) {
      progress.setTask("Error - check logs!");
      logger.log(Level.SEVERE, "Failed to load priors", e);
    } finally {
      if (reader != null) {
        try {
          reader.close();
        } catch (IOException e) { /* ignore */ }
      }
      if (conn != null) {
        try {
          conn.rollback();
        } catch (SQLException e) {/* ignore */}
        try {
          conn.close();
        } catch (SQLException e) {/* ignore */}
      }
      progress.complete();
    }
  }

  /**
   * Removes databases where the sequence files are missing.
   * Regenerates background files when they are missing.
   * @param ds the database source.
   * @param binDir the directory containing fasta-get-markov.
   * @param dbDir the directory containing the fasta databases.
   * @param status the status.
   */
  public static void purge(SQLiteDataSource ds, File binDir, File dbDir, boolean warnOnly, MultiSourceStatus status) throws InterruptedException {
    Connection conn = null;
    int regeneratedBgs = 0;
    int deletedFiles = 0;
    int deletedListings;
    int deletedCategories;
    Progress progress = null;
    try {
      progress = new Progress(status, "File check");
      // get connection to db
      conn = ds.getConnection();
      conn.setAutoCommit(false);
      // select fileSeq, fileBg from tblSequenceFile
      PreparedStatement pstmt = conn.prepareStatement(SQL.SELECT_ALL_SEQUENCE_FILES);
      PreparedStatement pstmtDelFile = conn.prepareStatement(SQL.DELETE_SEQUENCE_FILE);
      ResultSet rset = pstmt.executeQuery();
      while (rset.next()) {
        long fileID = rset.getLong(1);
        AlphStd alphabet = AlphStd.fromInt(rset.getInt(2));
        String seqFileName = rset.getString(3);
        String bgFileName = rset.getString(4);
        progress.setTask("Checking existence of " + seqFileName);
        // check existence of files in the dbDir with the names
        File seqFile = new File(dbDir, seqFileName);
        File bgFile = new File(dbDir, bgFileName);
        if (!warnOnly) {
          if (!seqFile.exists()) {
            // when sequence is missing remove the database entry and associated background
            progress.setTask("Removing entry for missing file " + seqFile.getName());
            logger.log(Level.FINE, "Deleting listing for missing sequence file " + seqFile);
            pstmtDelFile.setLong(1, fileID);
            deletedFiles += pstmtDelFile.executeUpdate();
            if (bgFile.exists() && !bgFile.delete()) {
              logger.log(Level.WARNING, "Could not delete unused file " + bgFile);
            }
          } else if (!bgFile.exists()) {
            // when background is missing run fasta-get-markov to regenerate it
            progress.setTask("Recalculating missing background for " + seqFile.getName(), 0, seqFile.length());
            logger.log(Level.FINE, "Regenerating missing background for sequence file " + seqFile);
            try {
              SequenceProcessor.generateBackground(binDir, alphabet, seqFile, bgFile, null);
              regeneratedBgs++;
            } catch (IOException e) {
              logger.log(Level.SEVERE, "Could not regenerate missing background file " + bgFile, e);
            }
          }
        } else {
          if (!seqFile.exists()) {
            logger.log(Level.WARNING, "Contains entry for missing sequence file " + seqFile);
          }
          if (!bgFile.exists()) {
            logger.log(Level.WARNING, "Contains entry for missing background file " + bgFile);
          }
        }
      }
      rset.close();
      pstmt.close();
      pstmtDelFile.close();
      // now check for and delete empty listings
      progress.setTask("Cleaning up unused listings");
      pstmt = conn.prepareStatement(SQL.DELETE_LISTING_WITHOUT_SEQUENCE_FILE);
      deletedListings = pstmt.executeUpdate();
      pstmt.close();
      // now check for and delete empty categories
      progress.setTask("Cleaning up unused categories");
      pstmt = conn.prepareStatement(SQL.DELETE_CATEGORY_WITHOUT_LISTING);
      deletedCategories = pstmt.executeUpdate();
      pstmt.close();
      // all done, commit changes
      progress.setTask("Commiting changes");
      conn.commit();
      conn.close();
      conn = null;
      logger.log(Level.INFO, "Purged database entries for " + deletedFiles +
          " missing files as well as " + deletedListings +
          " empty listings and " + deletedCategories + " empty categories. " +
          " Also regenerated " + regeneratedBgs + " background files.");
    } catch (SQLException e) {
      progress.setTask("Error - check logs!");
      logger.log(Level.SEVERE, "Failed to purge invalid database entries", e);
    } finally {
      if (conn != null) {
        try {
          conn.rollback();
        } catch (SQLException e) {/* ignore */}
        try {
          conn.close();
        } catch (SQLException e) {/* ignore */}
      }
      progress.complete();
    }
  }

  /**
   * Stops the downloader tasks by sending interrupts.
   */
  private static class ShutdownCmd implements CmdHandler {
    /** matches exit, stop or shutdown command */
    private static final Pattern EXIT_CMD_RE = Pattern.compile(
        "^(?:exit|stop|shutdown)(?:\\s.*)?$", Pattern.CASE_INSENSITIVE);
    /** Help message displayed to the user about the exit command */
    private static final CmdHelp EXIT_INFO =
        new CmdHelp("exit", "Stops downloads, finishes processing sequences and performs an orderly shutdown.");
    /** reference to the download threads so they can be interrupted */
    private ExecutorService downloadWorkers;

    /**
     * Constructor for the shutdown command.
     * @param downloadWorkers reference to the download threads so they can be interrupted.
     */
    public ShutdownCmd(ExecutorService downloadWorkers) {
      this.downloadWorkers = downloadWorkers;
    }

    /**
     * Test for shutdown commands (allows exit, stop or shutdown) and
     * activate the shutdown of the downloader tasks.
     * @param command the command that was read off the console
     * @return status mode if a handled command.
     */
    @Override
    public CmdResult command(String command) {
      Matcher m = EXIT_CMD_RE.matcher(command);
      if (m.matches()) {
        downloadWorkers.shutdownNow();
        return CmdResult.STATUS_MODE;
      } else {
        return CmdResult.UNHANDLED;
      }
    }

    /**
     * Gets the help for the shutdown command.
     * @return the help for the shutdown command.
     */
    @Override
    public List<CmdHelp> getHelp() {
      return Arrays.asList(EXIT_INFO);
    }
  }

  /**
   * Runs the passed updater classes
   * @param ds the database source, the database should be initialised already.
   * @param binDir the directory containing the fasta_get_markov program.
   * @param dbDir the directory containing the fasta databases.
   * @param updaters a list of classes to run to update the database.
   * @param status access to the command line to display status messages.
   */
  public static void update(SQLiteDataSource ds,
      File binDir, File dbDir, List<String> updaters, MultiSourceStatus status) {
    logger.log(Level.INFO, "Started update of sequence databases in " + dbDir);
    // check that we actually have some updaters selected
    if (updaters.size() > 0) {
      // this thread is responsible for unpacking files and generating the background
      // this might be fairly processor intensive so by limiting it to one thread
      // we stop the computer from overloading itself (hopefully)
      // possibly change the thread count to be configurable in future
      ExecutorService sequenceWorker = Executors.newSingleThreadExecutor();
      // these threads will download the sequences and queue them for processing
      // this should not be processor intensive so allow all to run in parallel
      // hopefully this won't cause bandwidth problems
      ExecutorService downloadWorkers = Executors.newFixedThreadPool(updaters.size());
      // mediate concurrent access to the database with a lock
      ReentrantReadWriteLock lock = new ReentrantReadWriteLock();
      // initialize and launch all the downloader threads
      List<Future<?>> downloaderTasks = new ArrayList<Future<?>>();
      for (String className : updaters) {
        try {
          // use some reflection to find the constructor
          Class<?> updaterClass = Class.forName(className);
          Constructor<?> updaterConstructor = updaterClass.getConstructor(
              SQLiteDataSource.class, ReentrantReadWriteLock.class, File.class,
              File.class, ExecutorService.class, MultiSourceStatus.class);
          Object updater = updaterConstructor.newInstance(ds, lock, binDir,
              dbDir, sequenceWorker, status);
          // launch the task
          downloaderTasks.add(downloadWorkers.submit((Callable<?>)updater));
        } catch (Exception e) {
          logger.log(Level.SEVERE, "Failed to start updater " + className, e);
        }
      }
      // we don't have any more updaters so allow threads to quit when they are done
      downloadWorkers.shutdown();
      // allow the user to shutdown the downloading threads
      status.registerCmdHandler("shutdown", new ShutdownCmd(downloadWorkers));
      // Wait for all the downloader tasks to complete.
      // Technically I could use awaitTermination but then I wouldn't be able to
      // access the exceptions from the tasks
      for (Future<?> future : downloaderTasks) {
        try {
          future.get(); // the result will be null on success so we don't check it
        } catch (InterruptedException e) {
          break;
        } catch (CancellationException e) { /* ignore */
        } catch (ExecutionException e) {
          logger.log(Level.WARNING, "Updater exception", e);
        }
      }
      // now that the workers are shutdown there is no reason to have this command registered.
      status.deregisterCmdHandler("shutdown");
      // just in case we got interrupted, make sure the workers exit
      downloadWorkers.shutdownNow();
      while (true) { // wait for all the downloader tasks to exit
        try {
          if (downloadWorkers.awaitTermination(1, TimeUnit.DAYS)) break;
        } catch (InterruptedException e) { /* ignore */ }
      }
      // At this point no threads should be adding sequences to be processed so
      // we can start the orderly shutdown of the worker
      sequenceWorker.shutdown(); // allow all sequences to be processed
      // wait for remaining sequences to be processed
      while (true) { // wait for the sequences to be processed
        try {
          if (sequenceWorker.awaitTermination(1, TimeUnit.DAYS)) break;
        } catch (InterruptedException e) { /* ignore */ }
      }
    }
    logger.log(Level.INFO, "Finished update of sequence databases in " + dbDir);
  }

  public static void markFilesObsolete(SQLiteDataSource ds, String globFileNames) {
    Connection conn = null;
    try {
      conn = ds.getConnection();
      PreparedStatement pstmt;
      pstmt = conn.prepareStatement(SQL.MARK_GLOB_SEQUENCE_FILES);
      pstmt.setInt(1, 1);
      pstmt.setString(2, globFileNames);
      int count = pstmt.executeUpdate();
      pstmt.close();
      conn.close(); conn = null;
      logger.log(Level.INFO, "Marked obsolete " + count + " sequence files matching glob: " + globFileNames);
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Failed to mark obsolete sequence files matching glob: " + globFileNames, e);
    } finally {
      if (conn != null) {
        try {
          conn.close();
        } catch (SQLException e) {/* ignore */}
      }

    }
  }

  /**
   * Removes all sequences marked as obsolete.
   * @param ds the database data source.
   * @param dbDir the directory containing the fasta files.
   */
  public static void removeObsolete(SQLiteDataSource ds, File dbDir) {
    Connection conn = null;
    try {
      conn = ds.getConnection();
      conn.setAutoCommit(false);
      PreparedStatement pstmt;
      ResultSet rset;
      pstmt = conn.prepareStatement(SQL.SELECT_ALL_OBSOLETE_SEQUENCE_FILES);
      rset = pstmt.executeQuery();
      while (rset.next()) {
        String seqFileName = rset.getString(1);
        String bgFileName = rset.getString(2);
        File seqFile = new File(dbDir, seqFileName);
        File bgFile = new File(dbDir, bgFileName);
        logger.log(Level.FINE, "Removing obsolete sequence file " + seqFile);
        if (seqFile.exists() && !seqFile.delete()) {
          logger.log(Level.WARNING, "Failed to delete obsolete sequence file " + seqFile);
        }
        logger.log(Level.FINE, "Removing obsolete sequence background file " + bgFile);
        if (bgFile.exists() && !bgFile.delete()) {
          logger.log(Level.WARNING, "Failed to delete obsolete sequence background file " + bgFile);
        }
      }
      rset.close();
      pstmt.close();
      pstmt = conn.prepareStatement(SQL.DELETE_OBSOLETE_SEQUENCE_FILES);
      int deletedFiles = pstmt.executeUpdate();
      pstmt.close();
      pstmt = conn.prepareStatement(SQL.DELETE_LISTING_WITHOUT_SEQUENCE_FILE);
      int deletedListings = pstmt.executeUpdate();
      pstmt.close();
      pstmt = conn.prepareStatement(SQL.DELETE_CATEGORY_WITHOUT_LISTING);
      int deletedCategories = pstmt.executeUpdate();
      pstmt.close();
      conn.commit();
      conn.close();
      conn = null;
      logger.log(Level.INFO, "Removed " + deletedFiles +
          " obsolete files as well as " + deletedListings +
          " empty listings and " + deletedCategories + " empty categories.");
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Failed to remove obsolete sequences", e);
    } finally {
      if (conn != null) {
        try {
          conn.rollback();
        } catch (SQLException e) {/* ignore */}
        try {
          conn.close();
        } catch (SQLException e) {/* ignore */}
      }
    }
  }

  /**
   * Replaces reserved values (comma and newline) in CSV with their HTML encodings.
   * @param value the string to be written as a CSV value.
   * @return the input string with instances of comma and newline replaced.
   */
  private static String replaceComma(String value) {
    return value.replace(",", "&#44;").replace("\n", "&#10;");
  }

  /**
   * Convert a boolean into a yes/no value.
   * @param value the boolean to convert.
   * @return the yes or no string respectively for true or false.
   */
  private static String y(boolean value) {
    return value ? "yes" : "no";
  }

  /**
   * Writes a record for all the sequence databases into a csv.
   * @param ds the database.
   * @param csv_file the target file to write the csv file.
   * @param status access to the command line to display status messages.
   */
  public static void generateCSV(SQLiteDataSource ds,
      File csv_file, MultiSourceStatus status) {
    logger.log(Level.INFO, "Generating CSV for sequence databases");
    Progress progress = new Progress(status, "Generating CSV");
    progress.setTask("Loading sequence database");
    boolean anyError = false;
    SequenceDBList sequenceDB = new SequenceDBList(ds);
    Writer out = null;
    try {
      progress.setTask("Writing CSV");
      out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(csv_file), "UTF-8"));
      out.append("#<db_root>,<prot?>,<nuc?>,<short_seqs?>,<db_menu_name>,<db_long_name>\n");
      try {
        progress.setTaskProgress(0, sequenceDB.countListings());
        for (Category category : sequenceDB.getCategories(false, EnumSet.allOf(AlphStd.class))) {
          logger.log(Level.FINE, "Writing category " + category.getName());
          out.append(",,,,-----").append(replaceComma(category.getName())).
              append("-----,,,\n");

          for (Listing listing : sequenceDB.getListings(category.getID())) {
            // list most recent version
            SequenceVersion version = sequenceDB.getVersions(listing.getID()).get(0);
            SequenceDB dnaSeqs = version.getSequenceFile(AlphStd.DNA);
            SequenceDB protSeqs = version.getSequenceFile(AlphStd.PROTEIN);
            String name = null;
            if (dnaSeqs != null) {
              name = dnaSeqs.getSequenceName();
            } else if (protSeqs != null) {
              name = protSeqs.getSequenceName();
            }
            if (name == null) continue;
            String baseName = name.substring(0, name.lastIndexOf('.'));
            boolean prot = protSeqs != null;
            boolean dna = dnaSeqs != null;
            boolean short_seqs = (dnaSeqs != null && dnaSeqs.getAverageLength() < 1000 && dnaSeqs.getMaxLength() < 5000);
            logger.log(Level.FINER, "Writing listing " + listing.getName());
            out.append(baseName).append(",").append(y(prot)).
                append(",").append(y(dna)).append(",").append(y(short_seqs)).
                append(",").append(replaceComma(listing.getName())).append(",").
                append(replaceComma(listing.getDescription())).append("\n");
            progress.addTaskPartDone(1);
          }
        }
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error accessing sequence database", e);
        progress.setTask("Error accessing sequence database - check logs!");
        anyError = true;
      }
      out.close(); out = null;
    } catch (IOException e) {
      logger.log(Level.SEVERE, "Error writing CSV file", e);
      progress.setTask("Error writing CSV file - check logs!");
      anyError = true;
    } finally {
      if (out != null) {
        try {
          out.close();
        } catch (IOException e) { /* ignore */ }
      }
    }

    if (anyError) {
      logger.log(Level.WARNING, "Errors occurred generating CSV for sequence databases");
    } else {
      logger.log(Level.INFO, "Successfully generated CSV for sequence databases");
    }
    progress.complete();
  }

  /**
   * Returns a binary string for a boolean.
   * @param value the boolean to convert into a string.
   * @return "1" for true and "0" for false.
   */
  private static String b(boolean value) {
    return value ? "1" : "0";
  }

  /**
   * Writes out a UTF-8 encoded string with a preceding length when given a Java string.
   * @param out the file to write to.
   * @param string the string to write.
   * @throws IOException if a error occurs writing the string.
   */
  private static void writeBinString(RandomAccessFile out, String string) throws IOException {
    byte[] utf8str = string.getBytes("UTF-8");
    out.writeInt(utf8str.length);
    out.write(utf8str);
  }

  /**
   * Write out a whole lot of strings in the binary format.
   * @param out the file to write the strings.
   * @param fields the strings to write.
   * @throws IOException if a error occurs writing a string.
   */
  private static void writeBinStrings(RandomAccessFile out, String ... fields) throws IOException {
    out.writeInt(fields.length);
    for (String field : fields) writeBinString(out, field);
  }

  /**
   * Write out a database listing in the binary format.
   * @param out the file to write the listing.
   * @param sequenceDB the sequence database.
   * @param listing the listing.
   * @throws SQLException if a query fails.
   * @throws IOException if writing out something fails.
   */
  private static void writeBinListing(RandomAccessFile out, SequenceDBList sequenceDB, Listing listing) throws SQLException, IOException {
    SequenceVersion version = sequenceDB.getVersions(listing.getID()).get(0);
    SequenceDB dnaSeqs = version.getSequenceFile(AlphStd.DNA);
    SequenceDB protSeqs = version.getSequenceFile(AlphStd.PROTEIN);
    String name = null;
    if (dnaSeqs != null) {
      name = dnaSeqs.getSequenceName();
    } else if (protSeqs != null) {
      name = protSeqs.getSequenceName();
    }
    if (name == null) return;
    String baseName = name.substring(0, name.lastIndexOf('.'));
    boolean prot = protSeqs != null;
    boolean dna = dnaSeqs != null;
    boolean short_seqs = (dnaSeqs != null && dnaSeqs.getAverageLength() < 1000 && dnaSeqs.getMaxLength() < 5000);
    writeBinStrings(out, baseName, b(prot), b(dna), b(short_seqs), listing.getName(), listing.getDescription());
  }

  /**
   * Write out a csv index (old sequence database format).
   * @param ds a database source.
   * @param index_file the file to write the csv index.
   * @param status status message handler.
   */
  public static void generateCSVIndex(SQLiteDataSource ds, File index_file, MultiSourceStatus status) {
    logger.log(Level.INFO, "Generating CSV index for sequence databases");
    Progress progress = new Progress(status, "Generating CSV index");
    progress.setTask("Loading sequence database");
    boolean anyError = false;

    SequenceDBList sequenceDB = new SequenceDBList(ds);
    RandomAccessFile out = null;
    try {
      final int INT_SIZE = 4;
      progress.setTask("Writing index");
      out = new RandomAccessFile(index_file, "rw");
      // truncate file
      out.setLength(0);
      // move to start
      out.seek(0);
      // write file format version
      out.writeInt(1);
      // write Category count and Listing count.
      int categoryCount = sequenceDB.countCategories();
      int totalListingCount = sequenceDB.countListings();
      progress.setTaskProgress(0, categoryCount + totalListingCount);
      out.writeInt(categoryCount);
      out.writeInt(totalListingCount);
      long categoryLoc = out.getFilePointer();
      long listingLoc = categoryLoc + categoryCount * INT_SIZE;
      // skip over category and listing index and write the category names
      out.seek((3 + categoryCount + totalListingCount) * INT_SIZE);
      // write the category names
      for (Category category : sequenceDB.getCategories(false, EnumSet.allOf(AlphStd.class))) {
        writeBinString(out, category.getName());
        progress.addTaskPartDone();
      }
      // write the category listings
      int entryIndex = 0;
      for (Category category : sequenceDB.getCategories(false, EnumSet.allOf(AlphStd.class))) {
        long loc = out.getFilePointer();
        if (loc > Integer.MAX_VALUE) throw new Error("Indexes larger than 4 GB can't be created due to use of 4 byte offsets");
        out.seek(categoryLoc);
        out.writeInt((int) loc);
        categoryLoc = out.getFilePointer();
        out.seek(loc);
        int listingCount = sequenceDB.countListings(category.getID());
        out.writeInt(entryIndex);
        out.writeInt(listingCount);
        writeBinString(out, category.getName());
        for (Listing listing : sequenceDB.getListings(category.getID())) {
          loc = out.getFilePointer();
          if (loc > Integer.MAX_VALUE)
            throw new Error("Indexes larger than 4 GB can't be created due to use of 4 byte offsets");
          out.seek(listingLoc);
          out.writeInt((int) loc);
          listingLoc = out.getFilePointer();
          out.seek(loc);
          writeBinListing(out, sequenceDB, listing);
          entryIndex++;
          progress.addTaskPartDone();
        }
      }
      out.close(); out = null;
    } catch (IOException e) {
      logger.log(Level.SEVERE, "Error writing CSV index file", e);
      progress.setTask("Error writing CSV index file - check logs!");
      anyError = true;
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Error accessing sequence database", e);
      progress.setTask("Error accessing sequence database - check logs!");
      anyError = true;
    } finally {
      if (out != null) {
        try {
          out.close();
        } catch (IOException e) { /* ignore */ }
      }
    }
    if (anyError) {
      logger.log(Level.WARNING, "Errors occurred generating CSV for sequence databases");
    } else {
      logger.log(Level.INFO, "Successfully generated CSV for sequence databases");
    }
    progress.complete();
  }

  /**
   * Convert a integer log level into the enum understood by the logger.
   * @param logLevel the log level as read from the command line argument.
   * @return the log level as understood by the logger.
   */
  private static Level getLevel(int logLevel) {
    switch (logLevel) {
      case 1: return Level.SEVERE;
      case 2: return Level.WARNING;
      case 3: return Level.INFO;
      case 4: return Level.CONFIG;
      case 5: return Level.FINE;
      case 6: return Level.FINER;
      case 7: return Level.FINEST;
    }
    return logLevel <= 0 ? Level.OFF : Level.ALL;
  }

  /**
   * Generate a file name for the log file.
   * @param dbDir the fasta database directory.
   * @return a file name for the log file.
   */
  private static String getLogPattern(File dbDir) {
    DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
    String today = dateFormat.format(new Date());
    // try to avoid naming conflicts when the updater is run multiple times in one day
    int run = 1;
    try {
      while (GlobFilter.find(dbDir, "update_" + today + "_" + run + ".*.log").length > 0) {
        run++;
        // do give up eventually
        if (run > 100) {
          run = 0;
          break;
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
      run = 0;
    }
    return new File(dbDir, "update_" + today + "_" + run).getPath() + ".%g.%u.log";
  }

  /**
   * Try to determine the bin directory by checking in order the existence of:
   * 1. the command line option --libexec
   * 2. the environment variable MEME_LIBEXC_DIR
   * 3. the configured value for the ${prefix}/bin
   * 4. the environment variable PATH looking for a folder containing the program fasta-get-markov
   * If options 1-3 are used then the folder is also checked for
   * fasta-get-markov and a IOException is thrown when it is missing.
   * @param config the parsed command line arguments.
   * @return the directory containing the program fasta-get-markov.
   * @throws IOException if fasta-get-markov cannot be found.
   */
  private static File getBinDir(JSAPResult config) throws IOException {
    File binDir = null;
    if (config.contains("libexec_dir")) {
      binDir = config.getFile("libexec_dir");
    } else if (System.getenv("MEME_LIBEXEC_DIR") != null) {
      binDir = new File(System.getenv("MEME_LIBEXEC_DIR"));
    } else {
      // look for MemeSuite.properties
      InputStream inStream = Thread.currentThread().getContextClassLoader().getResourceAsStream("MemeSuite.properties");
      if (inStream != null) {
        Properties properties = new Properties();
        properties.load(inStream);
        if (properties.containsKey("libexec.dir")) {
          binDir = new File(properties.getProperty("libexec.dir"));
        }
      }
    }
    found:
    if (binDir == null) {
      // search path
      String path_var = System.getenv("PATH");
      String[] exe_paths = path_var.split(File.pathSeparator);
      for (String exe_path : exe_paths) {
        File fastaGetMarkov = new File(exe_path, "fasta-get-markov");
        if (fastaGetMarkov.exists() && fastaGetMarkov.canExecute()) {
          binDir = new File(exe_path);
          break found;
        }
      }
      throw new FileNotFoundException("Cannot find runnable fasta-get-markov program in path.");
    } else {
      File fastaGetMarkov = new File(binDir, "fasta-get-markov");
      if (!fastaGetMarkov.exists() || !fastaGetMarkov.canExecute()) {
        throw new FileNotFoundException("Cannot find runnable fasta-get-markov program at location " + binDir);
      }
    }
    return binDir;
  }


  /**
   * Run the program.
   * @param args the program arguments.
   * @throws Exception when something goes wrong.
   */
  public static void main(String[] args) throws Exception {
    Properties setup = new Properties();
    setup.load(UpdateSequenceDB.class.getResourceAsStream("UpdateSequenceDB.properties"));

    List<String> ids = new ArrayList<String>();
    // now get a list of the downloader classes from the properties
    Pattern DL_CLASS_RE = Pattern.compile("^dl\\.([a-z0-9]+)\\.class$");
    for (String property : setup.stringPropertyNames()) {
      Matcher m = DL_CLASS_RE.matcher(property);
      if (m.matches()) ids.add(m.group(1));
    }
    Collections.sort(ids);

    List<String> updaters = new ArrayList<String>();
    List<Parameter> parameters = new ArrayList<Parameter>();

    for (String id : ids) {
      // add the class to the list of downloaders
      String qualifiedClass = setup.getProperty("dl." + id + ".class");
      updaters.add(qualifiedClass);
      // add an option to enable the downloader
      Switch aSwitch = new Switch(qualifiedClass).setLongFlag(id);
      String description = setup.getProperty("dl." + id + ".description", id);
      if (description != null) aSwitch.setHelp("Enable downloading " + description + ".");
      parameters.add(aSwitch);
      // add an option to disable the downloader
      aSwitch = new Switch("no_" + qualifiedClass).setLongFlag("no_" + id);
      description = setup.getProperty("dl." + id + ".description", id);
      if (description != null) aSwitch.setHelp("Disable downloading " + description + ".");
      parameters.add(aSwitch);
    }

    // define some special parsers for files
    StringParser FILE_PARSER = FileStringParser.getParser().setMustBeFile(true);
    StringParser DIR_PARSER = FileStringParser.getParser().setMustBeDirectory(true);

    // add the always existing switches
    parameters.add(
      new FlaggedOption("updater").setStringParser(new ClassnameStringParser())
          .setLongFlag("updater").setShortFlag('u')
          .setAllowMultipleDeclarations(true)
          .setHelp("Specify the classname of a custom updater [experimental].")
    );
    parameters.add(
        new FlaggedOption("obsolete")
            .setLongFlag("obsolete")
            .setAllowMultipleDeclarations(true)
            .setHelp("Mark all sequence databases matching the GLOB pattern file name(s) as obsolete and exit.")
    );
    parameters.add(
        new FlaggedOption("load_priors").setStringParser(FILE_PARSER)
            .setLongFlag("priors")
            .setHelp("Add priors to the database for the entries mentioned in the tab separated file and exit.")
    );
    parameters.add(
      new Switch("remove_obsolete")
          .setLongFlag("delete_old").setShortFlag('d')
          .setHelp("Sequence databases marked as obsolete (on a previous update) will be deleted.")
    );
    parameters.add(
        new Switch("retain_missing")
        .setLongFlag("retain_missing")
        .setHelp("Sequence databases that the file has been lost for will be kept in the database.")
    );
    parameters.add(
      new QualifiedSwitch("csv_dir").setStringParser(FILE_PARSER).
          setLongFlag("csv").setShortFlag('c').
          setHelp("Create a csv file and index file that lists all the databases to " +
              "enable backwards compatibility with older releases.")
    );
    parameters.add(
      new FlaggedOption("bin_dir").setStringParser(DIR_PARSER).
          setLongFlag("bin").setShortFlag('b').
          setHelp("Specify the path to the bin directory where the " +
              "fasta-get-markov tool can be located.")
    );
    parameters.add(
      new FlaggedOption("log_file").setStringParser(FILE_PARSER)
          .setLongFlag("log").setShortFlag('l')
          .setHelp("Specify the file that logging should be written to, " +
              "otherwise create a log file in db directory.")
    );
    parameters.add(
      new FlaggedOption("log_level").setStringParser(JSAP.INTEGER_PARSER)
          .setLongFlag("verbosity").setShortFlag('v').setDefault("3")
          .setHelp("Specify the logging level [1-8].")
    );
    parameters.add(
      new UnflaggedOption("db_dir").setStringParser(DIR_PARSER).setRequired(true)
          .setHelp("Specify the directory containing the databases.")
    );

    SimpleJSAP jsap = new SimpleJSAP(
      "update-sequence-db",
      "Update sequence databases; optionally remove old databases.",
      parameters.toArray(new Parameter[parameters.size()])
    );
    JSAPResult config = jsap.parse(args);
    if (jsap.messagePrinted()) {
      System.exit(1);
    }

    // See if any specified updaters are given.
    Boolean specified_only = false;
    Iterator<String> iter = updaters.iterator();
    while (iter.hasNext()) if (config.getBoolean(iter.next())) {
      specified_only = true; 
      break;
    }

    // Remove all but specified updaters if one or more are given.
    iter = updaters.iterator();
    if (specified_only) {
      while (iter.hasNext()) if (! config.getBoolean(iter.next())) iter.remove();
    } else {
      // remove any disabled updaters
      while (iter.hasNext()) if (config.getBoolean("no_" + iter.next())) iter.remove();
    }

    // add any additional updaters
    updaters.addAll(Arrays.asList(config.getStringArray("updater")));

    // determine bin dir
    File binDir;
    try {
      binDir = getBinDir(config);
    } catch (FileNotFoundException e) {
      System.err.println(e.getMessage());
      System.exit(1);
      return; // make compiler happy
    }

    File dbDir = config.getFile("db_dir");
    // create the database directory if it doesn't exist
    if (!dbDir.exists() && !dbDir.mkdirs()) {
      System.err.println("Unable to create database directory!");
      System.exit(1);
    }
    File log_dir = new File(dbDir, "logs");
    if (!log_dir.exists() && !log_dir.mkdir()) {
      System.err.println("Unable to create log directory!");
      System.exit(1);
    }

    // configure the logger
    Level level = getLevel(config.getInt("log_level"));
    String logPattern;
    if (config.contains("log_file")) {
      logPattern = config.getFile("log_file").getPath();
    } else {
      logPattern = getLogPattern(log_dir);
    }
    FileHandler handler = new FileHandler(logPattern);
    handler.setFormatter(new SimpleFormatter());
    logger.addHandler(handler);
    logger.setLevel(level);
    logger.setUseParentHandlers(false);

    // setup the interactive status
    MultiSourceStatus status = new MultiSourceStatus(250, 5000, 500);
    try {
      // load the database
      SQLiteDataSource ds = getInitialisedDatabase(dbDir);
      // are we marking obsolete databases?
      String[] obsoleteGlobs = config.getStringArray("obsolete");
      File priorsFile = config.getFile("load_priors");
      if (obsoleteGlobs != null && obsoleteGlobs.length > 0) {
        Progress progress = new Progress(status, "Marking files obsolete");
        for (String obsoleteGlob : obsoleteGlobs) {
          progress.setTask("GLOB " + obsoleteGlob);
          markFilesObsolete(ds, obsoleteGlob);
        }
        progress.complete();
      } else if (priorsFile != null) {
        loadPriors(ds, dbDir, priorsFile, status);
      } else {
        // purge missing databases
        purge(ds, binDir, dbDir, config.getBoolean("retain_missing", false), status);
        // remove obsolete databases
        if (config.getBoolean("remove_obsolete")) removeObsolete(ds, dbDir);
        // run the selected updaters
        update(ds, binDir, dbDir, updaters, status);
        // create the optional backwards compatibility views
        if (config.getBoolean("csv_dir")) {
          File csv_dir = config.getFile("csv_dir", dbDir);
          generateCSV(ds, new File(csv_dir, "fasta_db.csv"), status);
          generateCSVIndex(ds, new File(csv_dir, "fasta_db.index"), status);
        }
      }
      status.addStatusSource("Done").done();
    } finally {
      status.shutdownAndWait();
      System.out.println();
    }
  }

  /**
   * Parses a class name to ensure that the class exists.
   */
  private static class ClassnameStringParser extends StringParser {
    @Override
    public Object parse(String s) throws ParseException {
      try {
        Class.forName(s);
      } catch (ClassNotFoundException e) {
        throw new ParseException("Could not find custom updater class: " + s);
      }
      return null;  // generated code
    }
  }
}
