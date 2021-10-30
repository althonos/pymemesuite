#ifndef DISPLAY_H
#define DISPLAY_H

#include "meme.h"

extern void print_results(
  DATASET *dataset,      /* the dataset IN */
  MODEL *model,          /* the model */
  CANDIDATE *candidates, /* candidate models found IN */
  FILE* outfile          /* file for text output IN */
);

extern void record_results(
  DATASET *dataset,			         /* the dataset IN */
  MODEL *model,				           /* the model IN */
  MOTIF_SUMMARY *motif_summaries /* summaries of final motifs IN */
);

/**********************************************************************
  print_meme_file_xml

  Print MEME results in XML format. See DTD embeded in
  print_meme_header_xml for description of MEME document.
 **********************************************************************/
extern void print_meme_file_xml(
  MODEL *model,                   /* the model IN */
  DATASET *dataset,               /* the dataset IN */
  DATASET *neg_dataset,           // the control dataset IN
  int nmotifs,                    /* number of motifs IN */
  MOTIF_SUMMARY *motif_summaries, /* list of final motif properties IN */
  char *stopping_reason,          /* description of reason for stopping IN */
  char* xml_filename              /* full path to output file for xml IN */
);

/**********************************************************************
print_meme_file_html

Print MEME results in HTML format.
Format XML as HTML using a stylesheet and an XSLT
**********************************************************************/
extern void print_meme_file_html(
    char* stylesheet_file_path,   /* path to MEME XSL stylesheet IN */
    char* input_file_path,        /* path to XML input file IN */
    char* output_file_path        /* path to HTML output file IN */
);

extern void print_theta(
  int imotif,		 /* motif number */
  char *consensus,  // single-letter consensus of motif
  int format,		 /* 1 = floating point
              	    2 = integer */ 
  int nsites,    /* number of sites (discrete) */
  THETA theta,	 /* theta */
  int w,		     /* width of motif */
  double log_ev,	/* log motif E-value */
  char *str_space,	/* space for printing strand direction */
  DATASET *dataset,	/* the dataset */ 
  FILE *outfile	 	  /* output file */
);

extern void print_zij(
  DATASET *dataset,			/* the dataset */
  MODEL *model				/* the model */
);

extern void print_wij(
  DATASET *dataset		/* the dataset */
);

extern char *get_consensus(
  THETA theta,			/* motif theta */
  int w,			/* width of motif */
  DATASET *dataset,		/* the dataset */
  int N,			/* number of letters for each position */
  double min_prob		/* minimum cumulative prob for N letters */
);

extern void print_command_summary(
  MODEL *model,     /* the model IN */
  DATASET *dataset, /* the dataset IN */
  FILE *outfile     /* where to print IN */
);

extern void print_dataset_summary (
  DATASET *dataset, /* the dataset IN */
  FILE *outfile     /* where to print IN */
);

extern void print_summary(
  MODEL *model,     /* the model IN */
  DATASET *dataset, /* the dataset IN */
  int nmotifs,      /* number of motifs IN */
  FILE *outfile     /* where to print IN */
);

extern void print_meme_doc(
  FILE *outfile         /* where to print IN */
);

/*
  Get a single-letter consensus from a model.
*/
extern char *get_single_letter_consensus(
  MODEL *model,
  ALPH_T *alphabet
);

#endif

