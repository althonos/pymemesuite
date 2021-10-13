##
## Process this file with automake to produce Makefile.in
##

# Set installation location for doc files.
memedatadir = $(DESTDIR)${datadir}/${PACKAGE_NAME}-${PACKAGE_VERSION}
if WEBSITE
MEME_DOC_DIR = $(DESTDIR)$(WEBDIR)/doc
else
MEME_DOC_DIR = $(DESTDIR)$(memedatadir)/doc
endif
memelibdir = ${libdir}/${PACKAGE_NAME}-${PACKAGE_VERSION}
memelibexecdir = ${libexecdir}/${PACKAGE_NAME}-${PACKAGE_VERSION}


ARCHIVE_REVISION = $(shell cat $(top_srcdir)/ARCHIVE_REVISION)
ARCHIVE_DATE = $(shell cat $(top_srcdir)/ARCHIVE_DATE)

SEDSPEC = \
  -e 's%@ARCHIVE_REVISION@%$(ARCHIVE_REVISION)%' \
  -e 's%@ARCHIVE_DATE@%$(ARCHIVE_DATE)%' \
  -e 's%@ARCHIVE_VERSION@%$(VERSION)%' \
  -e 's%@TMP_DIR@%$(MEME_TEMP)%' \
  -e 's%@WHICHMONO@%$(MONO)%' \
  -e 's%@WHICHPERL@%$(PERL)%' \
  -e 's%@WHICHPYTHON@%$(PYTHON)%' \
  -e 's%@WHICHPYTHON3@%$(PYTHON3)%' \
  -e 's%@WHICHJAVA@%$(JAVA)%' \
  -e 's%@WHICHCONVERT@%$(CONVERT)%' \
  -e 's%@WHICHGHOSTSCRIPT@%$(GHOSTSCRIPT)%' \
  -e 's%@PERLLIBDIR@%$(memelibdir)/perl%' \
  -e 's%@PYTHONLIBDIR@%$(memelibdir)/python%' \
  -e 's%@JAVALIBDIR@%$(memelibdir)/java%' \
  -e 's%@LAM_BIN@%$(LAMBIN)%' \
  -e 's%@APPCONFIGDIR@%$(memedatadir)%' \
  -e 's%@BINDIR@%$(bindir)%' \
  -e 's%@LIBEXECDIR@%$(memelibexecdir)%' \
  -e 's%@MEMELOGS@%$(MEME_LOGS)%' \
  -e 's%@MEMEDB@%$(MEME_DB)%' \
  -e 's%@SITE_URL@%$(URL)%g' \
  -e 's%\@SOURCE_URL\@%$(SOURCE_URL)%g' \
  -e 's%@MAX_TIME@%$(MAXTIME)%' \
  -e 's%@DOCDIR@%$(MEME_DOC_DIR)/%' \
  -e 's%@maxmemory@%$(MAXMEMORY)%' \
  -e 's%@streme_length_factor@%$(STREME_LENGTH_FACTOR)%' \
  -e 's%@mpi_nprocs@%$(MPI_NPROCS)%'

BUILT_PERL_MOD = \
  Globals.pm \
  MemeWebUtils.pm \
  MotifUtils.pm \
  StatusPage.pm

NORMAL_PERL_MOD = \
  Alphabet.pm \
  Citation.pm \
  CheckingSAX.pm \
  DremeSAX.pm \
  ExecUtils.pm \
  HtmlMonolithWr.pm \
  JsonRdr.pm \
  JsonWr.pm \
  MastSAX.pm \
  MemeSAX.pm \
  MotifInDremeXML.pm \
  MotifInMemeXML.pm \
  MotifInStremeXML.pm \
  PriorUtils.pm \
  ReadFastaFile.pm \
  StremeSAX.pm \
  TomtomSAX.pm

OTHER_PERL_MOD = \
  CombDynProg.pm \
  DiffJSON.pm \
  DiffXML.pm \
  HypergeometricDynProg.pm

perlmoddir=$(memelibdir)/perl
perlmod_SCRIPTS = \
  $(BUILT_PERL_MOD) \
  $(NORMAL_PERL_MOD)

BUILT_PYTHON_LIB = \
  fasta-dinucleotide-shuffle.py \
  fasta-hamming-enrich.py

PYTHON_LIB = \
  alphabet.py \
  hypergeometric.py \
  sequence.py

PYTHON3_LIB = \
  alphabet_py3.py \
  hypergeometric_py3.py \
  sequence_py3.py

pythonlibdir=$(memelibdir)/python
pythonlib_SCRIPTS = \
  $(BUILT_PYTHON_LIB) \
  $(PYTHON_LIB) \
  $(BUILT_PYTHON3_LIB) \
  $(PYTHON3_LIB)

INPUT_FILES = $(wildcard *.in)

BUILT_UTILS = \
  ama-qvalues \
  beeml2meme \
  centrimo-plots \
  dreme_xml_to_html \
  dreme_xml_to_txt \
  elm2meme \
  fasta-center \
  fasta-dinucleotide-shuffle \
  fasta-fetch \
  fasta-hamming-enrich \
  fasta-make-index \
  fasta-most \
  fasta-re-match \
  fasta-subsample \
  fasta-unique-names \
  glam2html \
  glam2psfm \
  glam2scan2html \
  chen2meme \
  jaspar2meme \
  hart2meme-bkg \
  hartemink2psp \
  iupac2meme \
  mast_xml_to_html \
  mast_xml_to_txt \
  matrix2meme \
  meme-chip_html_to_tsv \
  meme-rename \
  meme_xml_to_html \
  nmica2meme \
  priority2meme \
  psp-gen \
  rsat-retrieve-seq \
  rsat-supported-organisms \
  rna2meme \
  scpd2meme \
  sites2meme \
  streme_xml_to_html \
  taipale2meme \
  tamo2meme \
  tomtom_xml_to_html \
  transfac2meme \
  uniprobe2meme \
  xstreme_html_to_tsv

PROGS = dreme \
  meme-chip \
  xstreme

BUILT_UTILS += \
  fasta-dinucleotide-shuffle \
  fasta-hamming-enrich

UTILS = \
  fasta-grep \
  plotgen \
  sd \
  $(BUILT_UTILS)

OBSOLETE_UTILS = \
  fasta-enriched-center \
  mhmm2html

# services
BUILT_DUALUSE_SERVICES = \
  tomtom_webservice
	
BUILT_WEBSERVICE_SERVICES = \
  ame_webservice \
  centrimo_webservice \
  tgene_webservice \
  dreme_webservice \
  fimo_webservice \
  glam2_webservice \
  glam2scan_webservice \
  gomo_webservice \
  mast_webservice \
  mcast_webservice \
  meme-chip_webservice \
  meme_webservice \
  momo_webservice \
  sea_webservice \
  spamo_webservice \
  streme_webservice \
  xstreme_webservice 

#utilities
WEBSITE_BUILT_WEBUTILS = \
  create-png-plots \
  cron-usage \
  lrun \
  make-usage-reports \
  texit \
  update-plot-usage \
  update-sequence-db

WEBSITE_WEBUTILS = \
  combine-meme-logs \
  plot-usage \
  plotgen \
  $(WEBSITE_BUILT_WEBUTILS)

MEME_SERVICES =
MEME_WEBUTILS =

# WEBSITE OR WEBSERVICE
if WEBSITE
MEME_SERVICES += $(BUILT_DUALUSE_SERVICES)
else !WEBSITE
if WEBSERVICE
MEME_SERVICES += $(BUILT_DUALUSE_SERVICES)
endif WEBSERVICE
endif !WEBSITE


if WEBSERVICE
MEME_SERVICES += $(BUILT_WEBSERVICE_SERVICES)
endif

if WEBSITE
MEME_WEBUTILS += $(WEBSITE_WEBUTILS)
endif


# Describe dependencies and build instructions for scripts.
# (substitute in local path to Perl).

# input files with unusual build names
MemeWebUtils.pm: MemeWebUtils.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@
MotifUtils.pm: MotifUtils.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@
StatusPage.pm: StatusPage.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@
Globals.pm: Globals.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@

# csh input files
%: %.csh.in
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@

# sh input files
%: %.sh.in
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@

# perl input files
%: %.pl.in
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@

if HAVE_PYTHON_3
dreme: dreme_py3.py.in alphabet_py3.py sequence_py3.py hypergeometric_py3.py ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
fasta-dinucleotide-shuffle: fasta-dinucleotide-shuffle-py3.in  ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
# Need to allow use as a Python library or for command line
# Python modules assume a file suffix of .py.
fasta-dinucleotide-shuffle.py: fasta-dinucleotide-shuffle
	cp fasta-dinucleotide-shuffle fasta-dinucleotide-shuffle.py
fasta-hamming-enrich: fasta-hamming-enrich-py3.in  ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
fasta-hamming-enrich.py: fasta-hamming-enrich
	cp  fasta-hamming-enrich fasta-hamming-enrich.py
glam2html: glam2html-py3.in ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
glam2psfm: glam2psfm-py3.in ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
glam2scan2html: glam2scan2html-py3.in ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
else
dreme: dreme.py.in alphabet.py sequence.py hypergeometric.py ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
# Need to allow use as a Python library or for command line
# Python modules assume a file suffix of .py.
fasta-dinucleotide-shuffle: fasta-dinucleotide-shuffle.py.in  ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
fasta-dinucleotide-shuffle.py: fasta-dinucleotide-shuffle
	cp fasta-dinucleotide-shuffle fasta-dinucleotide-shuffle.py
fasta-hamming-enrich: fasta-hamming-enrich.py.in  ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
fasta-hamming-enrich.py: fasta-hamming-enrich
	cp fasta-hamming-enrich fasta-hamming-enrich.py
glam2html: glam2html.py.in ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
glam2psfm: glam2psfm.py.in ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
glam2scan2html: glam2scan2html.py.in ../ARCHIVE_DATE ../ARCHIVE_REVISION
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
	rm -f ../src/$@
	ln $@ ../src/$@ # Copy to src for testing purposes
endif

# Set automake targets
bin_SCRIPTS = $(PROGS) 
memelibexec_SCRIPTS = $(UTILS) $(MEME_WEBUTILS) $(MEME_SERVICES)

# Define which files to delete if target is 'make clean'.
CLEANFILES = \
  *.pyc \
	__pycache__/* \
  $(BUILT_WEBSERVICE_SERVICES) \
  $(BUILT_PERL_MOD) \
  $(BUILT_PYTHON_LIB) \
  $(BUILT_PYTHON3_LIB) \
  $(BUILT_UTILS) \
  $(PROGS) \
  $(WEBSITE_BUILT_WEBUTILS) \
  $(MEME_SERVICES)

EXTRA_DIST = \
  $(NORMAL_PERL_MOD) \
  $(OTHER_PERL_MOD) \
  $(PYTHON_LIB) \
  $(PYTHON3_LIB) \
  combine-meme-logs \
  create-png-plots \
  dependencies.pl \
  fasta-grep \
  make-usage-reports \
  obo2dag.jar \
  plotgen \
  plot-usage \
  sd \
  $(INPUT_FILES)

