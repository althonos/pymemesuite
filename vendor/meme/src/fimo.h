#ifndef FIMO_H
#define FIMO_H

#include "alphabet.h"
#include "array-list.h"
#include "projrel.h"
#include "string-list.h"
#include "utils.h"

typedef enum {INVALID_THRESH, PV_THRESH, QV_THRESH} THRESHOLD_TYPE;

const char *threshold_type_to_string(THRESHOLD_TYPE type);

// Structure for tracking fimo command line parameters.
typedef struct options {

  bool allow_clobber;      // Allow overwritting of files in output directory.
  bool compute_qvalues;    // Compute q-values
  bool parse_genomic_coord;// Parse genomic coord. from seq. headers.
  bool text_only;          // Generate only plain text output
  bool max_strand;         // When scores available for both strands
                                // print only the max of the two strands.
  bool scan_both_strands;  // Scan forward and reverse strands
  bool skip_matched_sequence;  // Don't report matched sequence in --text mode

  char* bg_filename;            // Name of file file containg background freq.
  char* command_line;           // Full command line
  char* meme_filename;          // Name of file containg motifs.
  char* output_dirname;         // Name of the output directory
  char* seq_filename;           // Name of file containg sequences.
  char* seq_name;               // Use this sequence name in the output.

  int max_stored_scores; // Maximum number of matches to store per pattern.

  double alpha;       // Non-motif specific scale factor.
  double pseudocount; // Pseudocount added to Motif PSFM.
  double output_threshold; // Maximum p-value/q-value to report.

  ALPH_T *alphabet;    // Alphabet specified by MEME file.
  THRESHOLD_TYPE threshold_type;  // Type of output threshold.
  STRING_LIST_T* selected_motifs; // Indices of requested motifs.

  char *psp_filename; // Path to file containing position specific priors (PSP)
  char *prior_distribution_filename; // Path to file containing prior distribution
  char *pval_lookup_filename;   // Print p-value lookup table.
  char *html_path; // Path to FIMO HTML output file
  char *text_path; // Path to FIMO plain-text output file
  char *gff_path; // Path to FIMO GFF output file
  char *xml_path; // Path to FIMO XML output file
  char *cisml_path; // Path to CisML XML output file

  const char* HTML_FILENAME;    // Name of HTML output file.
  const char* TSV_FILENAME;     // Name of TSV output file.
  const char* GFF_FILENAME;     // Name of GFF output file.
  const char* XML_FILENAME;    // Name of FIMO XML output file.
  const char* CISML_FILENAME;   // Name of CisML XML output file.

  const char* usage; // Usage statment

} FIMO_OPTIONS_T;
#endif
