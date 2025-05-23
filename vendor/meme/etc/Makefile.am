##
## Process this file with automake to produce Makefile.in
##

AUTOMAKE_OPTIONS = no-dependencies

memelibexecdir = ${libexecdir}/${PACKAGE_NAME}-${PACKAGE_VERSION}

DATABASE_FILES = \
  db_general.csv \
  db_other_genomes.csv

SOURCE_CONFIG_FILES = \
  AME_config.xml.in \
  CENTRIMO_config.xml.in \
  DREME_config.xml.in \
  FIMO_config.xml.in \
  GLAM2_config.xml.in \
  GLAM2SCAN_config.xml.in \
  GOMO_config.xml.in \
  MAST_config.xml.in \
  MCAST_config.xml.in \
  MEMECHIP_config.xml.in \
  MEME_config.xml.in \
  MOMO_config.xml.in \
  SEA_config.xml.in \
  SPAMO_config.xml.in \
  STREME_config.xml.in \
  TGENE_config.xml.in \
  TOMTOM_config.xml.in \
  TOMTOM_SHORT_config.xml.in \
  XSTREME_config.xml.in

BUILT_CONFIG_FILES = \
  AME_$(VERSION).xml \
  CENTRIMO_$(VERSION).xml \
  DREME_$(VERSION).xml \
  FIMO_$(VERSION).xml \
  GLAM2_$(VERSION).xml \
  GLAM2SCAN_$(VERSION).xml \
  GOMO_$(VERSION).xml \
  MAST_$(VERSION).xml \
  MCAST_$(VERSION).xml \
  MEMECHIP_$(VERSION).xml \
  MEME_$(VERSION).xml \
  MOMO_$(VERSION).xml \
  SEA_$(VERSION).xml \
  SPAMO_$(VERSION).xml \
  STREME_$(VERSION).xml \
  TGENE_$(VERSION).xml \
  TOMTOM_$(VERSION).xml \
  TOMTOM_SHORT_$(VERSION).xml \
  XSTREME_$(VERSION).xml

SOURCE_ETC_FILES = \
  ame_template.html.in \
  beadstring.xsl.in \
  constants.js.in \
  constants.xsl.in \
  centrimo_template.html.in \
  dreme_template.html.in \
  gomo-to-html.xsl.in \
  logging.conf.in \
  mast_template.html.in \
  mcast_template.html.in \
  meme_template.html.in \
  meme-chip_template.html.in \
  momo_doc.js \
  momo_template.html.in \
  sea_template.html.in \
  spamo_template.html.in \
  streme_template.html.in \
  tgene_template.html.in \
  tomtom_template.html.in \
  xsl-wrapper.in \
  xstreme_template.html.in

NORMAL_ETC_FILES = \
  alphabet_utilities.js \
  ame_doc.js \
  ame_icon.png \
  ame_template.css \
  ame_template.js \
  block_diagram.css \
  block_diagram.js \
  centrimo_doc.js \
  centrimo_graph.js \
  centrimo_icon.png \
  centrimo_template.css \
  centrimo_template.js \
  citation.js \
  discovery_doc.js \
  cisml.css \
  cisml-to-gff.xsl \
  cisml-to-gff3.xsl \
  cisml-to-html.xsl \
  cisml-to-text.xsl \
  delay_draw.js \
  description_section.js \
  dreme_icon.png \
  dreme_template.css \
  dreme_template.js \
  eps_context.js \
  fimo_doc.js \
  gomo_doc.js \
  generate_js_wrapped_template.pl \
  generate_xsl_icons.pl \
  hmm-text.xsl \
  job_status.tmpl \
  mast_icon.png \
  mast_template.css \
  mast_template.js \
  mcast_icon.png \
  mcast_doc.js \
  mcast_template.css \
  mcast_template.js \
  meme.css \
  meme_icon.png \
  meme_template.css \
  meme_template.js \
  memechip_icon.png \
  meme-chip_doc.js \
  meme-chip_template.css \
  meme-chip_template.js \
  momo_icon.png \
  momo_template.css \
  momo_template.js \
  motif_logo.js \
  new_icon.png \
  prior1.plib \
  prior30.plib \
  sea_doc.js \
  sea_icon.png \
  sea_template.css \
  sea_template.js \
  set_maxtime.pl \
  simple-shared-doc.js \
  site_histogram.js \
  spamo_doc.js \
  spamo_graph.js \
  spamo_icon.png \
  spamo_template.css \
  spamo_template.js \
  streme_doc.js \
  streme_icon.png \
  streme_template.css \
  streme_template.js \
  submit_or_download_motif.js \
  template.eps \
  template.js \
  tgene_doc.js \
  tgene_icon.png \
  tgene_template.css \
  tgene_template.js \
  tomtom_doc.js \
  tomtom_icon.png \
  tomtom_template.css \
  tomtom_template.js \
  utilities.js \
  xstreme_doc.js \
  xstreme_icon.png \
  xstreme_template.css \
  xstreme_template.js

BUILT_ETC_FILES = \
  ame_template.html \
  beadstring.xsl \
  centrimo_template.html \
  citation.js.h \
  citation.js.pm \
  citation.js.xsl \
  constants.js \
  constants.xsl \
  constants.js.xsl \
  dreme_template.html \
  gomo-to-html.xsl \
  gomo_doc.js.xsl \
  icons.xsl \
  logging.conf \
  mast_template.html \
  mcast_template.html \
  meme.css.xsl \
  meme_template.html \
  meme-chip_template.html \
  momo_template.html \
  motif_logo_template.js \
  simple-shared-doc.js.xsl \
  sea_template.html \
  spamo_template.html \
  streme_template.html \
  tgene_template.html \
  tomtom_template.html \
  utilities.js.xsl \
  xstreme_template.html

%.h: %
	node $? C > $@

%.pm: %
	node $? perl > $@

%.xsl: %
	$(SED) 's%tname%$?%' $(srcdir)/xsl-wrapper.in | $(SED) '5r $?' > $@

SED_XML = -e 's%@SERVICE_VERSION@%$(VERSION)%' -e 's%@SERVICE_DIR@%$(memelibexecdir)%' -e 's%@drmaa_queue@%$(DRMAA_QUEUE)%' -e 's%@drmaa_queue_short@%$(DRMAA_QUEUE_SHORT)%' -e 's%@maxmemory@%$(MAXMEMORY)G%' -e 's%@maxmemory_short@%$(MAXMEMORY_SHORT)G%'

%_$(VERSION).xml: %_config.xml.in
	$(SED) $(SED_XML) $? | $(PERL) $(top_srcdir)/etc/set_maxtime.pl $(MAXTIME) $(MAXTIME_SHORT) > $@

SED_URLS = -e 's%\@SOURCE_URL\@%$(SOURCE_URL)%g' -e 's%\@SITE_URL\@%$(URL)%g' -e 's%\@DEV_CONTACT\@%$(DEV_CONTACT)%g' -e 's%\@AMIGO_URL\@%$(GOSERVER_URL)%'

constants.js: constants.js.in
	$(SED) $(SED_URLS) $? > $@

constants.xsl: constants.xsl.in
	$(SED) $(SED_URLS) $? > $@

beadstring.xsl: beadstring.xsl.in
	$(SED) $(SED_URLS) $? > $@

%.html: %.html.in
	$(SED) $(SED_URLS) $< > $@

%-to-html.xsl: %-to-html.xsl.in
	$(SED) $(SED_URLS) $< > $@

icons.xsl: ../doc/images/gomo_icon.png
	$(PERL) $(top_srcdir)/etc/generate_xsl_icons.pl $? > $@

motif_logo_template.js: template.eps
	$(PERL) $(top_srcdir)/etc/generate_js_wrapped_template.pl motif_logo_template $? > $@

logging.conf: logging.conf.in 
	$(SED) $(SED_URLS) $? > $@

memedatadir = $(datadir)/$(PACKAGE_NAME)-${PACKAGE_VERSION}
memedata_DATA = $(DATABASE_FILES) $(NORMAL_ETC_FILES) $(BUILT_ETC_FILES)

if WEBSERVICE
deploydir=$(OPAL_DEPLOY)

deploy_DATA = $(BUILT_CONFIG_FILES)
endif
 
CLEANFILES = $(BUILT_CONFIG_FILES) $(BUILT_ETC_FILES)

EXTRA_DIST = $(DATABASE_FILES) $(SOURCE_ETC_FILES) $(NORMAL_ETC_FILES) $(SOURCE_CONFIG_FILES) citation.js.pm citation.js.h
