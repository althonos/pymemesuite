#####################################################################
## FILE: Makefile.am
# AUTHOR: William Stafford Noble, Charles E. Grant
# CREATE DATE: 2003
# PROJECT: Meta-MEME
# COPYRIGHT: 2003, UW
#####################################################################

## Automake uses this file to generate a Makefile

#####################################################################
#  This is the Makefile for the documentation directory.
#####################################################################
SUBDIRS = examples images css js

copyright.html: copyright.html.in ../COPYING Makefile
	$(SED) -e "/INSERT_COPYING_FILE/r ../COPYING" -e "/INSERT_COPYING_FILE/d" $< > copyright.html

SOURCE_HTML_FILES = \
  download.html.in \
  download-index.html.in \
  general-faq.html.in \
  install.html.in \
  scripting.html.in 

BUILT_HTML_FILES = \
  download.html \
  download-index.html \
  general-faq.html \
  install.html \
  scripting.html

# Files Included in distribution and installed.
NORMAL_FILES = \
  ../etc/prior30.plib \
  alphabets.html \
  alphabet-format.html \
  alphtype.html \
  ama.html \
  ama-qvalues.html \
  ame.html \
  ame-output-format.html \
  ame-tutorial.html \
  authors.html \
  beadstring.html \
  bed-format.html \
  beeml2meme.html \
  bfile-format.html \
  centrimo.html \
  centrimo-output-format.html \
  centrimo-tutorial.html \
  tgene.html \
  tgene-tutorial.html \
  ceqlogo.html \
  chen2meme.html \
  tgene-output-format.html \
  cite.html \
  clustalw2fasta.html \
  clustalw2phylip.html \
  clustalw-format.html \
  compute-prior-dist.html \
  compute-uniform-priors.html \
  copyright.html \
  create-priors.html \
  dmix-format.html \
  draw-mhmm.html \
  dreme.html \
  dreme-tutorial.html \
  dust.html \
  elm2meme.html \
  fasta-center.html \
  fasta-dinucleotide-shuffle.html \
  fasta-fetch.html \
  fasta-format.html \
  fasta-get-markov.html \
  fasta-grep.html \
  fasta-hamming-enrich.html \
  fasta-holdout-set.html \
  fasta-io.html \
  fasta-make-index.html \
  fasta-most.html \
  fasta-shuffle-letters.html \
  fasta-subsample.html \
  fasta-unique-names.html \
  fimo.html \
  fimo-output-format.html \
  fimo-tutorial.html \
  fisher_exact.html \
  fitevd.html \
  gendb.html \
  getsize.html \
  glam2.html \
  glam2format.html \
  glam2mask.html \
  glam2scan.html \
  glam2_alphabet.html \
  glam2_tut.html \
  GLAM2_method.pdf \
  godag-format.html \
  gomo.html \
  gomo-output-format.html \
  gomo_highlight.html \
  ismb94.pdf \
  iupac.html \
  iupac2meme.html \
  jaspar2meme.html \
  log-hmm.html \
  mast.html \
  matrix2meme.html \
  mcast.html \
  mcast-output-format.html \
  meme.html \
  meme2alph.html \
  meme-chip.html \
  meme-chip-output-format.html \
  meme-format.html \
  meme-get-motif.html \
  meme-rename.html \
  meme2images.html \
  meme2meme.html \
  mhmme.html \
  mhmm.html \
  mhmm2html.html \
  mhmm-format.html \
  mhmmscan.html \
  mhmmscan-format.html \
  mhmms.html \
  momo.html \
  momo-output-format.html \
  motif-shuffle-columns.html \
  motif-consensus.html \
  motif_conversion.html \
  motiph.html \
  news.txt \
  nmica2meme.html \
  overview.html \
  pmp_bf.html \
  priority2meme.html \
  psm-format.html \
  psp-format.html \
  psp-gen.html \
  purge.html \
  qvalue.html \
  reconcile-tree-alignment.html \
  reduce-alignment.html \
  release-notes.html \
  remove-alignment-gaps.html \
  rna2meme.html \
  scpd2meme.html \
  server_down.html \
  shadow.html \
  sea.html \
  sea-output-format.html \
  sites2meme.html \
  spamo.html \
  spamo-output-format.html \
  streme.html \
  streme-output-format.html \
  streme-tutorial.html \
  taipale2meme.html \
  tamo2meme.html \
  tomtom.html \
  tomtom-output-format.html \
  transfac2meme.html \
  transfac-format.html \
  uniprobe2meme.html \
  update-sequence-db.html \
  xstreme.html

SED_URLS = -e 's%\@SOURCE_URL\@%$(SOURCE_URL)%g' -e 's%\@SITE_URL\@%$(URL)%g' -e 's%\@DEV_CONTACT\@%$(DEV_CONTACT)%g'

%.html: %.html.in
	$(SED) $(SED_URLS) $< > $@

memedocdir = ${datadir}/${PACKAGE_NAME}-${PACKAGE_VERSION}/doc
export memedocdir
memedoc_DATA =  $(NORMAL_FILES) $(BUILT_HTML_FILES)

# force the generation of copyright.html when making a website
BUILT_SOURCES = copyright.html

CLEANFILES = $(BUILT_HTML_FILES) copyright.html

EXTRA_DIST = copyright.html.in $(NORMAL_FILES) $(SOURCE_HTML_FILES)
