##
## Process this file with automake to produce Makefile.in
##

##############################################

AUTOMAKE_OPTIONS = no-dependencies


SUBDIRS = \
  compute_prior_dist_example_output_files \
  sample_opal_scripts

SEDSPEC = \
  -e 's%@WHICHPERL@%$(PERL)%'

EXAMPLES = \
  glam2-aa-alph.txt \
  glam2-dna-alph.txt \
  sample-dna-aln.glam2 \
  sample-protein-aln.glam2 \
  sample-dna-motif.meme \
  sample-protein-motif.meme \
  small-tide-mod.txt

BUILT_PERL = \
  create_examples

# perl input files
create_examples: create_examples.pl.in
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@

noinst_SCRIPTS = $(BUILT_PERL)

# load the EXAMPLE_OUTPUT_FILES variable
include examples.mk

# load the EXAMPLE_DATASETS variable
include example-datasets.mk

examplesdir = ${memedocdir}/examples
export examplesdir
nobase_dist_examples_DATA = \
	$(EXAMPLES) \
	$(EXAMPLE_DATASETS) \
	$(EXAMPLE_OUTPUT_FILES) \
	examples.mk \
	example-datasets.mk \
	create_examples.pl.in 

CLEANFILES = create_examples 
