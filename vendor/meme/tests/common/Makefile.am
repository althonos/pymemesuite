## $Id: Makefile.am 1421 2006-11-01 01:15:40Z cegrant $
##
## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = no-dependencies

SUBDIRS = motif.crp0 motif.Puf3p

DATAFILES = \
	adh.s \
	alphabet_dna.txt \
	alphabet_rna.txt \
	At.s \
	chr22.fasta \
	crp0.meme.html \
	crp0.meme.xml \
	crp0.meme.txt \
	crp0.fasta \
	crp0.s \
	crp0-test.fasta \
	dna_with_covalent.alph \
	dna_with_covalent.bg \
	dna_with_covalent.fa \
	dna_with_covalent.meme \
	dros.fasta \
	dros.oops.html \
	extended_dna.alph \
	farntrans5.s \
	INO_up800.s \
	Jaspar-subset.meme \
	Jaspar-subset2.meme \
	Klf1-200-100.fa \
	Klf1-200-100.combined.fa \
	Klf1-200-100.combined.nofasta.fa \
	Klf1-200-100.no_fasta.fa \
	Klf1-200-100-shuffled.py2.fa \
	Klf1-200-100-shuffled.py3.fa \
	Klf1-200.fa \
	Klf1.fa \
	lex0.s \
	lipocalin.s \
	lipo.fasta \
	lipo.meme.html \
	lipo.meme.old.html \
	lipo-test.fasta \
	MCM1.fasta \
	mini-drosoph.s \
	non-uniform.bg \
	noq.s \
        nrmix.100bp.fa \
	nr.freq \
	nt.freq \
	psp_test.psp \
	psp_test.s \
	Puf3p-200.s \
	Puf3p-20.s \
	Puf3p.s \
	sample.meme \
	simple.aln \
	simple.meme \
	test.aln \
	test_tree.txt \
	uprot.freq \
	yeast.nc.1.freq \
	yeast.nc.3.freq \
	yeast.nc.6.freq 

if HAVE_PYTHON_3
Klf1-200-100-shuffled.fa: Klf1-200-100-shuffled.py3.fa
	cp Klf1-200-100-shuffled.py3.fa $(top_builddir)/tests/common/Klf1-200-100-shuffled.fa
else
Klf1-200-100-shuffled.fa: Klf1-200-100-shuffled.py2.fa
	cp Klf1-200-100-shuffled.py2.fa $(top_builddir)/tests/common/Klf1-200-100-shuffled.fa
endif

check_SCRIPTS = Klf1-200-100-shuffled.fa

EXTRA_DIST = $(DATAFILES)
CLEANFILES = Klf1-200-100-shuffled.fa
