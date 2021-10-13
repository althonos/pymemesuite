package au.edu.uq.imb.memesuite.updatedb;

import au.edu.uq.imb.memesuite.util.GlobFilter;
import au.edu.uq.imb.memesuite.db.*;
import com.martiansoftware.jsap.*;
import com.martiansoftware.jsap.stringparsers.FileStringParser;
import org.sqlite.SQLiteConfig;
import org.sqlite.SQLiteDataSource;

import java.io.*;
import java.lang.Runtime;
import java.sql.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Date;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import static au.edu.uq.imb.memesuite.db.SQL.*;

/**
 * Updates selected sequence databases with FASTA index files
 */
public class FastaIndexer {

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
   * Adds the fileSeqIndex column to the tblSequenceFile table 
   * if it doesn't already exist.
   * @param ds the DataSource.
   * @throws SQLException when the SQL doesn't work.
   */
  private static void updateSequenceFileTable(SQLiteDataSource ds) throws SQLException {
    logger.log(Level.INFO, "Update tbleSequenceFile columns.");
    Connection connection = null;
    try {
      connection = ds.getConnection();
      Statement statement = connection.createStatement();
      logger.log(Level.INFO, "Adding fileSeqIndex column to tblSequenceFile");
      ResultSet rs = statement.executeQuery(TEST_INDEX_FILE_COLUMN);
      String schema = rs.getString(1);
      if (!schema.contains("fileSeqIndex")) {
        statement.executeUpdate(ADD_FILE_SEQ_INDEX_COLUMN);
      }
      statement.close();
      connection.close();
      connection = null;
    } 
    catch (SQLException e) {
      logger.log(Level.WARNING, "Failed to create fileSeqIndex column. " + e);
    } 
    finally {
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) {
          logger.log(Level.WARNING, "Failed to close connection to database.");
        }
      }
    }
    logger.log(Level.INFO, "Successfully updated tblSequenceFile table.");
  }

  /**
   * For each genome database create an index for the FASTA file.
   * @param dbDir the name of the database directory.
   * @param ds the DataSource.
   * @throws SQLException when the SQL doesn't work.
   */
  private static void indexFastaFiles(String dbDir, SQLiteDataSource ds) throws SQLException {
    logger.log(Level.INFO, "Update tbleSequenceFile columns.");
    Connection connection = null;
    try {
      connection = ds.getConnection();
      Statement select = connection.createStatement();
      PreparedStatement update = connection.prepareStatement(UPDATE_FASTA_INDEX_FILE);

      ResultSet rs = select.executeQuery(SELECT_ENSEMBL_DNA_FASTA_FILES);
      while (rs.next()) {
        Integer id = rs.getInt(1);
        String fastaFile = rs.getString(2);
        String indexFile = rs.getString(3);
        if (indexFile == null) {
          String cmd = "index-fasta-file " + dbDir + "/" + fastaFile + " " 
            + dbDir + "/" + fastaFile + ".fai";
          System.out.println("Creating index for " + dbDir + fastaFile);
          Process ps = null;
          try {
            ps = Runtime.getRuntime().exec(cmd);
            ps.waitFor();
          }
          catch (Exception e) {
            logger.log(Level.WARNING, "Failed to create index for " + dbDir + fastaFile + "." );
            logger.log(Level.WARNING,  e.toString());
          }
          if (ps.exitValue() != 0) {
            logger.log(Level.WARNING, "Failed to create index for " + dbDir + fastaFile + ".");
            BufferedReader errorReader = new BufferedReader(
            new InputStreamReader(ps.getErrorStream()));
            String line;
            try {
              while ((line = errorReader.readLine()) != null) {
                System.out.println(line);
              }
              errorReader.close();
            }
            catch (Exception e) {
              logger.log(Level.WARNING, e.toString());
            }
          }
          else {
            update.setString(1, fastaFile + ".fai");
            update.setInt(2, id);
            int numChanged = update.executeUpdate();
            if (numChanged != 1) {
              logger.log(Level.WARNING, "Error updating record " + id.toString() + "." );
            }
          }
        }
      }
      select.close();
      update.close();
      connection.close();
      connection = null;
    } 
    catch (SQLException e) {
      logger.log(Level.WARNING, "Failed to create fileSeqIndex column. " + e);
    } 
    finally {
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) {
          logger.log(Level.WARNING, "Failed to close connection to database.");
        }
      }
    }
    logger.log(Level.INFO, "Successfully updated tblSequenceFile table.");
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
    SQLiteDataSource ds = getDBSource(db);
    updateSequenceFileTable(ds);
    return ds;
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
    setup.load(UpdateSequenceDB.class.getResourceAsStream("FastaIndexer.properties"));

    List<Parameter> parameters = new ArrayList<Parameter>();

    // define some special parsers for files
    StringParser FILE_PARSER = FileStringParser.getParser().setMustBeFile(true);
    StringParser DIR_PARSER = FileStringParser.getParser().setMustBeDirectory(true);

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
      "fasta-indexer",
      "Generate index files for genomic databases.",
      parameters.toArray(new Parameter[parameters.size()])
    );
    JSAPResult config = jsap.parse(args);
    if (jsap.messagePrinted()) {
      System.exit(1);
    }

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
    if (!dbDir.exists()) {
      System.err.println("Database directory doesn't exit!");
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

    SQLiteDataSource ds = getInitialisedDatabase(dbDir);
    indexFastaFiles(dbDir.toString(), ds);
  }
}
