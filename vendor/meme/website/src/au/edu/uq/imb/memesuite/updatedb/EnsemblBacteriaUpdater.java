package au.edu.uq.imb.memesuite.updatedb;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.util.MultiSourceStatus;
import org.apache.commons.net.PrintCommandListener;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.sqlite.SQLiteDataSource;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This updates the Ensembl Bacteria Genomes databases
 */
public class EnsemblBacteriaUpdater extends SequenceUpdater {
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.updatedb.ensembl.bacteria");
  private static Pattern dnaFilePattern1 = Pattern.compile("^.*\\.dna_sm\\.primary_assembly\\.fa\\.gz$");
  private static Pattern dnaFilePattern2 = Pattern.compile("^.*\\.dna_sm\\.toplevel\\.fa\\.gz$");
  private static Pattern aaFilePattern = Pattern.compile("^.*\\.pep\\.all\\.fa\\.gz$");
  private static Pattern releasePattern = Pattern.compile("/release-(\\d+)/");
  private static final String GENOMES = "Ensembl Bacteria Genomes and Proteins";
  private static final int ENSEMBL_RETRIEVER_ID = 3;
  private static final int RETRY_COUNT = 3;
  private static final int ERROR_COUNT = 10;
  private String scrape;
  private String host;
  private String ftpLog;

  public EnsemblBacteriaUpdater(SQLiteDataSource dataSource,
      ReentrantReadWriteLock dbLock, File binDir, File dbDir,
      ExecutorService worker, MultiSourceStatus statusManager) {
    super("Ensembl Bacteria Genomes Downloader", dataSource, dbLock, binDir, dbDir, worker, statusManager);
    Properties conf = loadConf(EnsemblBacteriaUpdater.class, dbDir, "EnsemblBacteriaUpdater.properties");
    scrape = conf.getProperty("scrape", "http://bacteria.ensembl.org/info/website/ftp/index.html").trim();
    host = conf.getProperty("host", "ftp.ensemblgenomes.org").trim();
    ftpLog = conf.getProperty("ftp.log", "").trim();
  }

  protected List<EnsemblGenome> scrapeAvailableGenomes() {
    progress.setTask("Querying available genomes", 0, -1);
    List<EnsemblGenome> genomes = new ArrayList<EnsemblGenome>();
    // get a listing of species
    try {
      Document doc = Jsoup.connect(scrape).get();
      Elements rows = doc.select("table.data_table tbody tr");
      logger.log(Level.INFO, "Found " + rows.size() + " rows in genome table.");
      for (Element row : rows) {
        Elements rowItems = row.select("td");
        String scientificName = rowItems.get(0).select("a").get(0).text();
        String dnaPath = new URL(rowItems.get(1).select("a").attr("href")).getPath();
        String proteinPath = new URL(rowItems.get(3).select("a").attr("href")).getPath();
        logger.log(Level.FINE, "Adding Ensembl DNA and protein genomes for " + scientificName);
        genomes.add(new EnsemblGenome(scientificName, scientificName, AlphStd.DNA, false, dnaPath, dnaFilePattern1, dnaFilePattern2));
        genomes.add(new EnsemblGenome(scientificName, scientificName, AlphStd.PROTEIN, false, proteinPath, aaFilePattern));
      }
    } catch (IOException e) {
      logger.log(Level.SEVERE, "Unable to get species listing");
    }
    return genomes;
  }

  protected boolean determineFtpSource(FTPClient ftp, EnsemblGenome genome) throws IOException {
    progress.setTask("Determining full URL for " + genome, 0, -1);
    PatternFileFilter filter = new PatternFileFilter(genome.getFilePattern());
    logger.log(Level.INFO, "Looking for " + genome + " sequences at " + genome.getRemoteDir());
    FTPFile sequenceFile = filter.best(ftp.listFiles(genome.getRemoteDir(), filter));
    if (sequenceFile == null) {
      logger.log(Level.WARNING, "Skipping " + genome + " as no sequence found.");
      return false;
    }
    logger.log(Level.INFO, "Using file " + sequenceFile.getName() + " for " + genome);
    genome.setRemoteInfo(sequenceFile.getName(), sequenceFile.getSize());
    return true;
  }

  protected boolean determineFtpSource(FTPClient ftp, EnsemblGenome genome, int retryCount) throws IOException, InterruptedException {
    for (int attempt = 1; attempt <= retryCount; attempt++) {
      try {
        // this could happen if we're retrying
        if (!ftp.isConnected()) {
          // try to recreate the connection
          if (!loginFtp(ftp, host)) throw new IOException("Unable to log in to " + host);
        }
        return determineFtpSource(ftp, genome);
      } catch (IOException e) {
        logger.log(Level.WARNING, "Failed to determine file for " + genome + " " + attempt + " of " + retryCount, e);
        if (attempt == retryCount) throw e;
        if (ftp.isConnected()) {
          try {
            ftp.disconnect();
          } catch (IOException e2) { /* ignore */ }
        }
        Thread.sleep(TimeUnit.SECONDS.toMillis(10));
      }
    }
    return false;
  }

  @Override
  public Void call() {
    int errorCount = 0;
    FTPClient ftp = null;
    try {
      logger.log(Level.INFO, "Starting Ensembl Bacteria Genomes update");
      List<EnsemblGenome> genomes = scrapeAvailableGenomes();
      if (genomes.isEmpty()) {
        logger.log(Level.SEVERE, "Failed parsing the genome releases from the FAQ... No genomes found!");
        return null;
      }
      progress.setTask("Excluding Pre-exisiting Genomes", 0, -1);
      // query the database to see which genomes we already have
      // and remove them so we don't download them again
      Iterator<EnsemblGenome> iterator = genomes.iterator();
      while (iterator.hasNext()) {
        EnsemblGenome genome = iterator.next();
        if (sourceExists(genome, true)) {
          iterator.remove();
          logger.log(Level.INFO, "Already have " + genome);
        }
      }
      if (genomes.isEmpty()) return null;

      ftp = new FTPClient();
      if (!ftpLog.isEmpty()) {
        ftp.addProtocolCommandListener(new PrintCommandListener(new PrintWriter(new FileOutputStream(ftpLog)), true));
      }
      if (!loginFtp(ftp, host)) return null;
      for (EnsemblGenome genome : genomes) {
        try {
          checkWorkerTasks();
          if (!determineFtpSource(ftp, genome, RETRY_COUNT)) continue;
          checkWorkerTasks();
          if (!downloadFtpSource(ftp, genome, true, RETRY_COUNT)) continue;
          enqueueSequences(new EnsemblSequenceProcessor(genome));
        } catch (IOException e) {
          logger.log(Level.WARNING, "Skipped " + genome + " due to ftp errors", e);
          errorCount++;
          if (errorCount >= ERROR_COUNT) throw new IOException("Too many IO Exceptions", e);
        }
      }
      ftp.logout();
      ftp.disconnect();
      progress.complete();
      waitForWorkerTasks();
      logger.log(Level.INFO, "Finished Ensembl Bacteria Genomes update");
    } catch (ExecutionException e) { // only thrown by sequence processor
      cancelWorkerTasks();
      logger.log(Level.SEVERE, "Abandoning Ensembl Bacteria update due to failure to process sequences!", e);
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Abandoning Ensembl Bacteria update!", e);
    } catch (IOException e) {
      logger.log(Level.SEVERE, "Abandoning Ensembl Bacteria update!", e);
    } catch (InterruptedException e) {
      logger.log(Level.WARNING, "Ensembl Bacteria update interrupted!");
    } catch (RuntimeException e) {
      logger.log(Level.SEVERE, "RuntimeException!", e);
      throw e;
    } catch (Error e) {
      logger.log(Level.SEVERE, "Error!", e);
      throw e;
    } finally {
      if (ftp != null && ftp.isConnected()) {
        try {
          ftp.logout();
        } catch (IOException e) { /* ignore */ }
        try {
          ftp.disconnect();
        } catch (IOException e) { /* ignore */ }
      }
      progress.complete();
    }
    return null;
  }

  private class EnsemblSequenceProcessor extends SequenceProcessor {
    private EnsemblGenome genome;

    public EnsemblSequenceProcessor(EnsemblGenome genome) {
      super(dataSource, dbLock, binDir, dbDir, status);
      this.genome = genome;
    }

    @Override
    public void process() throws IOException, SQLException, InterruptedException {
      if (!unpackSequences(genome)) return;
      processSequences(genome);
      recordSequences(genome);
    }
  }

  private class EnsemblGenome extends AbstractFtpSource {
    private String commonName;
    private String scientificName;
    private AlphStd alphabet;
    private boolean abinitio;
    private String ftpDir;
    private String ftpName;
    private long ftpSize;
    private Pattern[] filePattern;
    // derived from ftpDir
    private String version;
    private int edition;
    // files
    private File packedFile;
    private File sequenceFile;
    private File bgFile;
    // stats
    private long seqCount;
    private long seqMinLen;
    private long seqMaxLen;
    private double seqAvgLen;
    private double seqStdDLen;
    private long seqTotalLen;

    public EnsemblGenome(String commonName, String scientificName,
        AlphStd alphabet, boolean abinitio, String path, Pattern... filePattern) {
      Matcher m = releasePattern.matcher(path);
      if (!m.find()) throw new IllegalArgumentException("Path does not include release version!");
      this.commonName = commonName;
      this.scientificName = scientificName;
      this.alphabet = alphabet;
      this.abinitio = abinitio;
      this.ftpDir = path;
      this.filePattern = filePattern;
      version = m.group(1);
      edition = Integer.parseInt(version, 10);
    }

    public String toString() {
      return "ENSMBL v" + version + " " + commonName + " (" +
          (abinitio ? "Ab Initio Predicted " : "") + alphabet + ")";
    }

    public Pattern[] getFilePattern() {
      return filePattern;
    }

    @Override
    public String getRemoteHost() {
      return host;
    }

    @Override
    public String getRemoteDir() {
      return ftpDir;
    }

    public void setRemoteInfo(String ftpName, long ftpSize) {
      this.ftpName = ftpName;
      this.ftpSize = ftpSize;
    }

    @Override
    public String getRemoteExt() {
      String ext;
      if (ftpName.endsWith(".fa.gz")) {
        ext = ".fa.gz";
      } else if (ftpName.endsWith(".tar.gz")) {
        ext = ".tar.gz";
      } else if (ftpName.endsWith(".zip")) {
        ext = ".zip";
      } else {
        ext = super.getRemoteExt();
      }
      return "." + alphabet.toString().toLowerCase() + ext;
    }

    @Override
    public String getRemoteName() {
      return ftpName;
    }

    @Override
    public long getRemoteSize() {
      return ftpSize;
    }

    @Override
    public int getRetrieverId() {
      return ENSEMBL_RETRIEVER_ID;
    }

    @Override
    public String getCategoryName() {
      return GENOMES;
    }

    @Override
    public String getListingName() {
      return commonName;
    }

    @Override
    public AlphStd guessAlphabet() {
      return alphabet;
    }

    @Override
    public long getSequenceEdition() {
      return edition;
    }

    @Override
    public String getSequenceVersion() {
      return version;
    }

    @Override
    public String getSequenceDescription() {
      String readmeURL = "ftp://" + getRemoteHost() + getRemoteDir() + "README";
      return "Downloaded from <a href=\"" + getRemoteUrl() + "\">" +
          getRemoteUrl() + "</a>. Refer to the dataset <a href=\"" +
          readmeURL + "\">README</a> for further details.";
    }

    @Override
    public String getNamePrefix() {
      String name = scientificName;
      name = name.replaceAll("[^a-zA-Z0-9\\p{Z}\\s_\\.-]", "");
      name = name.replaceAll("[\\p{Z}\\s]+", "_");
      return "ensembl_" + name + "_" + version;
    }

    @Override
    public String getListingDescription() {
      return "Ensembl Bacteria Genomes - " + commonName + " (<i>" + scientificName + "</i>) - " + version;
    }


    public void setSourceFile(File file) {
      packedFile = file;
    }

    @Override
    public List<File> getSourceFiles() {
      return Collections.singletonList(packedFile);
    }

    @Override
    public void setSequenceFile(File file) {
      sequenceFile = file;
    }

    @Override
    public File getSequenceFile() {
      return sequenceFile;
    }

    @Override
    public void setBgFile(File file) {
      bgFile = file;
    }

    @Override
    public File getBgFile() {
      return bgFile;
    }

    @Override
    public void setStats(long count, long minLen, long maxLen, double avgLen,
        double stdDLen, long totalLen) {
      this.seqCount = count;
      this.seqMinLen = minLen;
      this.seqMaxLen = maxLen;
      this.seqAvgLen = avgLen;
      this.seqStdDLen = stdDLen;
      this.seqTotalLen = totalLen;
    }

    @Override
    public long getSequenceCount() {
      return seqCount;
    }

    @Override
    public long getTotalLength() {
      return seqTotalLen;
    }

    @Override
    public long getMinLength() {
      return seqMinLen;
    }

    @Override
    public long getMaxLength() {
      return seqMaxLen;
    }

    @Override
    public double getAverageLength() {
      return seqAvgLen;
    }

    @Override
    public double getStandardDeviationLength() {
      return seqStdDLen;
    }
  }
}
