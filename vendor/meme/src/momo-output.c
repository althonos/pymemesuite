#include "motif-in.h"
#include "config.h"
#include "momo-output.h"
#include "io.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "ceqlogo.h"
#include "momo-algorithm.h"
#include "momo-modl.h"
#include "fisher_exact.h"
#include "html-monolith.h"

/**********************************************************************
 * This function saves MOMO results as a MEME motif file.
 *********************************************************************/
static void print_momo_text_file(FILE *text_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  fprintf(text_file, "MEME version %s\n\n", VERSION);
  fprintf(text_file, "Alphabet= %.*s\n\n", 20, summary.alph->symbols + 1);
  fprintf(text_file, "Background letter frequencies\n");
  
  const char* alph_letters = summary.alph_letters;
  ARRAY_T * bg_freqs = summary.bg_freqs;
  int i, j, k, l;
  
  for (i = 0; i < strlen(alph_letters); ++i) {
    fprintf(text_file, "%c %8.6f ", alph_letters[i], get_array_item(i, bg_freqs));
  }
  fprintf(text_file, "\n\n");

  ARRAYLST_T * mod_table_keys = summary.mod_table_keys;
  
  for (i = 0; i < arraylst_size(mod_table_keys); i++) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    ARRAYLST_T* motifs = mod_entry->motifinfos;
    for (j = 0; j < arraylst_size(motifs); ++j) {
      MOTIF_INFO_T* motifinfo = arraylst_get(j, motifs);
      MOTIF_T* motif = motifinfo->motif;

      fprintf(text_file, "MOTIF %s    %s\n", get_motif_id(motif), get_motif_id2(motif));
      unsigned long motif_list_size = (unsigned long) arraylst_size(motifinfo->seqs);
      double m1, e1;
      exp10_logx(motif->log_evalue/log(10.0), m1, e1, 1);
      fprintf(text_file, "letter-probability matrix: alength= %d w= %d nsites= %lu E= %3.1fe%+04.0f\n", (int) (strlen(alph_letters)), options.width, motif_list_size, m1, e1);

      MATRIX_T* freqs = get_motif_freqs(motif);
      for (k = 0; k < options.width; k++) {
	for (l = 0; l < strlen(alph_letters); l++) {
	  fprintf(text_file, "%8.6f\t", get_matrix_cell(k, l, freqs));
	}
	fprintf(text_file, "\n");
      }
      fprintf(text_file, "\n");
    }
  }
};

void momo_finish_tsv_file(FILE *tsv_file, MOMO_OPTIONS_T options) {
  // Finish the TSV output
  char *version_message = "# MoMo (Modification Motifs): Version " VERSION " compiled on " __DATE__ " at " __TIME__ "\n";
  fprintf(tsv_file, "\n%s", version_message);
  fprintf(tsv_file, "# The format of this file is described at %s/%s.\n", SITE_URL, "doc/momo-output-format.html");
  fprintf(tsv_file, "# %s\n", options.command_line);
  fclose(tsv_file);
}

void print_algorithm_name(MOMO_OPTIONS_T options, FILE *fp) {
  ALGORITHM_T algorithm = options.algorithm;
  if (algorithm == Simple) {
    fprintf(fp, "simple");
  } else if (algorithm == Motifx) {
    fprintf(fp, "motif-x");
  } else if (algorithm == Modl) { 
    fprintf(fp, "MoDL");
  } else {
    fprintf(fp, "UNKNOWN");
  }
}

void momo_print_motifs(
  JSONWR_T *json_output,
  FILE *tsv_file, 
  MOMO_OPTIONS_T options, 
  SUMMARY_T summary
) {
  ARRAYLST_T * mod_table_keys = summary.mod_table_keys;
  int i, j, k, m;
  bool pvalues_accurate = (!options.db_background && options.algorithm == Motifx);
  const char* alph_letters = summary.alph_letters;

  // Loop over mods.
  jsonwr_bool_prop(json_output, "pvalues_accurate", pvalues_accurate);
  jsonwr_property(json_output, "mods");
  jsonwr_start_array_value(json_output);
  for (i = 0; i < arraylst_size(mod_table_keys); ++i) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    char *mod_name = mod_entry->mod_name;
    ARRAYLST_T * motifs = mod_entry->motifinfos;

    // Skip mods with no motifs.
    if (arraylst_size(motifs) == 0) continue;

    // Loop over motifs.
    jsonwr_start_object_value(json_output);
    jsonwr_str_prop(json_output, "mod_name", mod_name);
    jsonwr_property(json_output, "motifs");
    jsonwr_start_array_value(json_output);
    for (j = 0; j < arraylst_size(motifs); ++j) {
      MOTIF_INFO_T* currmotifinfo = arraylst_get(j, motifs);
      MOTIF_T* currmotif = currmotifinfo->motif;
      char* motifid = currmotif->id + 1;

      char* motif_name = mm_malloc(strlen(motifid) + 2);
      char *motif_regexp = NULL;
      // Add a central 'X' in the motif name if doing single motif per mass.
      for (k=m=0; k < strlen(motifid); k++) {
        // Check for single motif per mass (may be "__" or "_number").
        if ( (k == 0) && (motifid[0]-'0' >= 0 && motifid[0]-'0' <= 9) ) {
          // Add missing central 'X' to name.
          motif_name[m++] = 'X';
          motif_name[m++] = motifid[0];
	} else if (motifid[k] == '_' && (motifid[k+1] == '_' || (motifid[k+1]-'0' >= 0 && motifid[k+1]-'0' <= 9))) { 
          // Convert '_' to 'X' for (missing) central residue.
          motif_name[m++] = '_';
          motif_name[m++] = 'X';
	} else {
	  motif_name[m++] = motifid[k];
	}
      }
      motif_name[m] = '\0';
      strncpy(currmotif->id+1, motif_name, MAX_MOTIF_ID_LENGTH);

      if (options.algorithm == Motifx || options.algorithm == Modl) {
        motif_regexp = mm_malloc(strlen(motifid) + 2);
        // Convert motif_name to PERL regular expression.
	for (k=m=0; k < strlen(motif_name); k++) {
	  if (motif_name[k] == 'x' || motif_name[k] == 'X') {
	    motif_regexp[m++] = '.';
	  } else if (motif_name[k]-'A' >= 0 && motif_name[k]-'Z' <= 0) {
	    motif_regexp[m++] = motif_name[k];
	  } else if (motif_name[k] == '[' || motif_name[k] == ']') {
	    motif_regexp[m++] = motif_name[k];
	  }
	}
	motif_regexp[m] = '\0';
      } else {				// algorithm = simple
        // Create expanded motif name.
        int w = options.width;
        int idw = strlen(currmotif->id);
        motif_name = mm_realloc(motif_name, idw + w + 4);
        for (k=0; k < w/2; k++) {
          motif_name[k] = 'x';
        }
        motif_name[k++] = '_';
        for (m=1; m < idw; m++, k++) {
          motif_name[k] = currmotif->id[m];
        }
        motif_name[k++] = '_';
        for (m=0; m < w/2; m++, k++) {
          motif_name[k] = 'x';
        }
        motif_name[k] = '\0';

        // Get regular expression for motif
        // space for maximal size regexp [all letters-1]...[all letters-1] X width + Null
        int alen = strlen(alph_letters);
        int re_size = options.width*(alen+1) + 1;
        motif_regexp = mm_malloc(re_size);
        motif_regexp[0] = '\0';
        char *letters = mm_malloc(alen + 1);
        MATRIX_T* freqs = get_motif_freqs(currmotif);
        int p = 0;
	for (k = 0; k < options.width; k++) {
          int n = 0;
	  for (m = 0; m < alen; m++) {
            double f = get_matrix_cell(k, m, freqs);
	    if (f != 0) letters[n++] = alph_letters[m];
	  }
          letters[n] = '\0';
          if (n == alen) {
            motif_regexp[p++] = '.';
          } else if (n == 1) {
            motif_regexp[p++] = letters[0];
          } else {
            motif_regexp[p++] = '[';
            strncpy(motif_regexp+p, letters, re_size-p);
            p += n;
            motif_regexp[p++] = ']';
          }
          motif_regexp[p] = '\0';
	}
        free(letters);
      }
      strncpy(currmotif->id+1, motif_name, MAX_MOTIF_ID_LENGTH);
      strncpy(currmotif->id2, motif_regexp, MAX_MOTIF_ID_LENGTH);

      char* motif_png = mm_malloc(strlen(motif_name) + 5);
      motif_png[0] = '\0';
      strncat(motif_png, motif_name, strlen(motif_name));
      strncat(motif_png, ".png", 4);
      // Replace [] wight bd because square brackets are significant in URI.
      int n; 
      for (n=0; n<strlen(motif_name); n++) {
        if (motif_png[n] == '[') { 
          motif_png[n] = 'b';
        } else if (motif_png[n] == ']') {
          motif_png[n] = 'd';
        }
      }
      
      char* logo_file = mm_malloc(strlen(options.output_dirname) + strlen(motif_name) + 2);
      logo_file[0] = '\0';
      strncat(logo_file, options.output_dirname, strlen(options.output_dirname));
      strncat(logo_file, "/", 1);
      strncat(logo_file, motif_png, strlen(motif_name));
      CL_create1(currmotif, false, false, "MoMo", logo_file, false, true);
      
      double score = currmotifinfo->score;
      double n_tests = currmotifinfo->n_tests;
      int fg_matches=0, fg_size=0, bg_matches=0, bg_size=0;
      double fold=0;
      double m1=0, m2=0, m3=0;
      double e1=0, e2=0, e3=0;
      double log_pvalue=0, log_norm_pvalue=0, log_evalue=0;

      if (options.algorithm == Motifx || options.algorithm == Modl) {
        // To emulate original motif-x we use the counts in the remaining sequences;
	// by default we use all the matches in all the sequences.
        if (options.algorithm == Motifx && options.harvard) {
          fg_matches = currmotifinfo->fg_matches;
          fg_size = currmotifinfo->fg_size;
          bg_matches = currmotifinfo->bg_matches;
          bg_size = currmotifinfo->bg_size;
        } else {
	  fg_matches = currmotifinfo->afg_matches;
	  fg_size = currmotifinfo->afg_size;
	  bg_matches = currmotifinfo->abg_matches;
	  bg_size = currmotifinfo->abg_size;
        }
        fold = (bg_matches*fg_size != 0) ? (double) fg_matches*bg_size / (bg_matches*fg_size) : 0;

        // Get p-value, adjusted p-value and E-value.

        // Calculate significance if motif-x.
	log_pvalue = getLogFETPvalue(fg_matches, fg_size, bg_matches, bg_size, false);
	exp10_logx(log_pvalue/log(10.0), m1, e1, 1);
        if (pvalues_accurate || options.printp) {
          double log_n_tests = log(n_tests);
	  log_norm_pvalue = LOGEV(log_n_tests, log_pvalue);
	  exp10_logx(log_norm_pvalue/log(10.0), m2, e2, 1);
	  log_evalue = log_pvalue + log_n_tests;
	  exp10_logx(log_evalue/log(10.0), m3, e3, 1);
          if (!pvalues_accurate) {
	    log_norm_pvalue = log_evalue = 0;
          }
        }
        currmotif->log_evalue = log_evalue;

        // Print TSV header if first motif.
	if (i==0 && j==0) fprintf(tsv_file, 
          "mod\tmotif\tregexp\tscore\tfg_match\tfg_size\tbg_match\tbg_size\tfg/bg"
          "\tunadjusted_p-value%s",
          (pvalues_accurate ? "\ttests\tadjusted_p-value\n" : "\n"));
        // Print TSV.
	fprintf(tsv_file, 
          "%s\t%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%.1f\t%3.1fe%+04.0f",
          mod_name, motif_name, motif_regexp, score, fg_matches, fg_size, bg_matches, bg_size, fold, m1, e1
        );
        if (pvalues_accurate) { 
	  fprintf(tsv_file, "\t%.0f\t%3.1fe%+04.0f\n", n_tests, m2, e2);
        } else {
	  fprintf(tsv_file, "\n");
        }

        // Print the p-values for debugging.
	if (options.printp && options.algorithm == Motifx) { 
          printf("%s score: %.2f ntests: %.0f p-value: %4.2fe%+04.0f norm_p-value: %4.2fe%+04.0f E-value: %4.2fe%+04.0f p-value_accurate? %s\n", 
            motif_name, score, n_tests, m1, e1, m2, e2, m3, e3, pvalues_accurate ? "yes" : "no");
        }
      } else if (options.algorithm == Simple) {
        currmotif->log_evalue = 0;
        fg_matches = currmotifinfo->fg_matches;
        // Print TSV header if first motif.
	if (i==0 && j==0) fprintf(tsv_file, "mod\tmotif\tregexp\tfg_matches\n");
        // Print TSV.
        fprintf(tsv_file, "%s\t%s\t%s\t%d\n", mod_name, motif_name, motif_regexp, fg_matches);
      }
      
      // Output JSON for motif
      jsonwr_start_object_value(json_output);
	jsonwr_str_prop(json_output, "motif_name", motif_name);
	jsonwr_str_prop(json_output, "motif_png", motif_png);
	jsonwr_str_prop(json_output, "motif_regexp", motif_regexp);
	jsonwr_dbl_prop(json_output, "score", score);
	jsonwr_lng_prop(json_output, "fg_matches", fg_matches);
	jsonwr_lng_prop(json_output, "fg_size", fg_size);
	jsonwr_lng_prop(json_output, "bg_matches", bg_matches);
	jsonwr_lng_prop(json_output, "bg_size", bg_size);
	jsonwr_dbl_prop(json_output, "fold", fold);
	jsonwr_dbl_prop(json_output, "m1", m1);
	jsonwr_dbl_prop(json_output, "e1", e1);
	jsonwr_lng_prop(json_output, "n_tests", (int) n_tests);
	jsonwr_dbl_prop(json_output, "m2", m2);
	jsonwr_dbl_prop(json_output, "e2", e2);
	jsonwr_property(json_output, "occurrences");
          jsonwr_start_array_value(json_output);
	    for (k = 0; k < arraylst_size(currmotifinfo->seqs); ++k) {
	      char *curr_motifinfo_seq = arraylst_get(k, currmotifinfo->seqs);
	      jsonwr_str_value(json_output, arraylst_get(k, currmotifinfo->seqs));
	    }
          jsonwr_end_array_value(json_output);
      jsonwr_end_object_value(json_output);

      // cleanup
      myfree(motif_name);
      myfree(motif_regexp);
      myfree(motif_png);
      myfree(logo_file);
    } // motif
    jsonwr_end_array_value(json_output);

    // Print out MoDL log
    if (options.algorithm == Modl) {

      MATRIX_T* bg_freqs = NULL;
      bg_freqs = get_count_matrix(bg_freqs, mod_entry->bg_seq_list, NULL, &options, &summary);
      for (j = 0; j < options.width; ++j) {
        for (k = 0; k < strlen(summary.alph_letters); ++k) {
          set_matrix_cell(j,k,get_matrix_cell(j,k,bg_freqs)/arraylst_size(mod_entry->bg_seq_list),bg_freqs);
        }
      }
      
      ARRAYLST_T* modl_ops = mod_entry->modl_ops;
      ARRAYLST_T* temp_list = arraylst_create();
      double minDL = INFINITY;
      double minstep = 0;
      double initDL = 0;
      jsonwr_property(json_output, "modl_log");
      jsonwr_start_object_value(json_output);
      jsonwr_property(json_output, "modl_steps");
      jsonwr_start_array_value(json_output);
      for (j = 0; modl_ops && (j < arraylst_size(modl_ops)); ++j) {
        jsonwr_start_object_value(json_output);
        MODL_STEP_T* step = arraylst_get(j, modl_ops);
	jsonwr_lng_prop(json_output, "step", j);
	jsonwr_dbl_prop(json_output, "score", step->score);
        do_step(step, temp_list, bg_freqs, options.max_motifs, &options, &summary, mod_entry);
        if (j == 0 ) initDL = step->score;
        if (step->score < minDL) {
          minDL = step->score;
          minstep = j;
        }
        jsonwr_property(json_output, "reg_exps");
        jsonwr_start_array_value(json_output);
        for (k = 0; k < arraylst_size(temp_list); ++k) {
	  jsonwr_str_value(json_output, regexmotif_to_string(arraylst_get(k, temp_list), mod_entry, &summary, &options));
        }
        jsonwr_end_array_value(json_output);
        jsonwr_end_object_value(json_output);
      }
      jsonwr_end_array_value(json_output);
      jsonwr_dbl_prop(json_output, "final_step", minstep);
      jsonwr_dbl_prop(json_output, "final_dl", minDL);
      jsonwr_dbl_prop(json_output, "decrease", minDL != INFINITY ? initDL-minDL : 0);
      jsonwr_end_object_value(json_output);

      free_matrix(bg_freqs);
    } // MoDL 
    jsonwr_end_object_value(json_output);

  } // mods
 jsonwr_end_array_value(json_output);

} // momo_print_motifs

/**********************************************************************
 * This function saves MOMO results as an HTML file and a TSV file
 *********************************************************************/
void print_momo_html_and_tsv_files(
  int argc, 
  char **argv,
  FILE *tsv_file,
  MOMO_OPTIONS_T options, 
  SUMMARY_T summary
) {
  int i = 0;
  HTMLWR_T *html_output = NULL;
  JSONWR_T *json_output = NULL;

  // setup the html output writer
  if ((html_output = htmlwr_create(get_meme_data_dir(), options.TEMPLATE_FILENAME, false))) {
    htmlwr_set_dest_name(html_output, options.output_dirname, options.HTML_FILENAME);
    htmlwr_replace(html_output, "momo_data.js", "data");
    json_output = htmlwr_output(html_output);
    if (json_output == NULL) die("Template does not contain data section.\n");
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "Failed to open HTML template file.\n");
    html_output = NULL;
    json_output = NULL;
  }

  // write out some information
  jsonwr_str_prop(json_output, "version", VERSION);
  jsonwr_str_prop(json_output, "revision", REVISION);
  jsonwr_str_prop(json_output, "release", ARCHIVE_DATE);
  jsonwr_str_prop(json_output, "program", "MoMo");
  jsonwr_args_prop(json_output, "cmd", argc, argv);

  // options
  jsonwr_property(json_output, "options");
  jsonwr_start_object_value(json_output);
  jsonwr_str_prop(json_output, "algorithm", algorithm_names[(int) options.algorithm]);
  if (options.psm_type) jsonwr_str_prop(json_output, "psm_type", options.psm_type);
  if (options.sequence_column) jsonwr_str_prop(json_output, "sequence_column", options.sequence_column);
  jsonwr_str_prop(json_output, "filetype", filetype_names[options.filetype]);
  jsonwr_lng_prop(json_output, "width", options.width);
  jsonwr_lng_prop(json_output, "seed", options.seed);
  jsonwr_bool_prop(json_output, "db_background", options.db_background);
  if (options.protein_database_filename) {
    jsonwr_str_prop(json_output, "protein_database", options.protein_database_filename);
    jsonwr_str_prop(json_output, "protein_database_format", filetype_names[options.bg_filetype]);
  }
  if (options.filter_field) {
    jsonwr_property(json_output, "filter");
    jsonwr_start_object_value(json_output);
    jsonwr_str_prop(json_output, "filter_field", options.filter_field);
    jsonwr_str_prop(json_output, "filter_type", filtertype_names[(int) options.filter_type]);
    jsonwr_dbl_prop(json_output, "filter_threshold", options.filter_threshold);
    jsonwr_end_object_value(json_output);
  }
  jsonwr_bool_prop(json_output, "remove_unknowns", options.remove_unknowns);
  jsonwr_lng_prop(json_output, "eliminate_repeat_width", options.eliminate_repeat_width);
  jsonwr_lng_prop(json_output, "min_occurrences", options.min_occurrences);
  jsonwr_bool_prop(json_output, "single_motif_per_mass", options.single_motif_per_mass);
  jsonwr_lng_prop(json_output, "hash_fasta_width", options.hash_fasta_width);
  if (options.algorithm == Motifx) {
    jsonwr_dbl_prop(json_output, "score_threshold", options.score_threshold);
    jsonwr_bool_prop(json_output, "harvard", options.harvard);
  } else if (options.algorithm == Modl) {
    jsonwr_lng_prop(json_output, "max_motifs", options.max_motifs);
    jsonwr_lng_prop(json_output, "max_iterations", options.max_iterations);
    jsonwr_lng_prop(json_output, "max_no_decrease", options.max_no_decrease);
  }
  jsonwr_end_object_value(json_output);

  // PTM files
  jsonwr_property(json_output, "ptm_files");
  jsonwr_start_array_value(json_output);
    ARRAYLST_T* ptm_filenames = options.phospho_filenames;
    for (i = 0; i < arraylst_size(ptm_filenames); ++i) {
      char* ptm_filename = arraylst_get(i, ptm_filenames);
      //jsonwr_start_object_value(json_output);
      //jsonwr_str_prop(json_output, "ptm_filename", ptm_filename); 
      //jsonwr_end_object_value(json_output);
      jsonwr_str_value(json_output, ptm_filename); 
    }
  jsonwr_end_array_value(json_output);

  // Summary data
  jsonwr_property(json_output, "summary");
  jsonwr_start_object_value(json_output);
  jsonwr_lng_prop(json_output, "num_mod", summary.num_mod); 
  jsonwr_lng_prop(json_output, "num_modtype", summary.num_modtype); 
  jsonwr_lng_prop(json_output, "num_mod_passing", summary.num_mod_passing); 
  jsonwr_lng_prop(json_output, "num_modtype_passing", summary.num_modtype_passing); 
  jsonwr_lng_prop(json_output, "num_bg_mod", summary.num_bg_mod); 
  jsonwr_end_object_value(json_output);

  // The motifs
  momo_print_motifs(json_output, tsv_file, options, summary);

  // Finish the HTML output.
  if (html_output) {
    if (htmlwr_output(html_output) != NULL) {
      die("Found another JSON replacement!\n");
    }
    htmlwr_destroy(html_output);
    html_output = NULL;
  }
  momo_finish_tsv_file(tsv_file, options);

} // print_momo_html_and_tsv_files

static void create_directory(MOMO_OPTIONS_T options) {
  
  const bool PRINT_WARNINGS = false;
  if (create_output_directory(
                              options.output_dirname,
                              options.allow_clobber,
                              PRINT_WARNINGS
                              )
      ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }
}

/**********************************************************************
 * This function saves the MOMO results as a set of files in a
 * directory:
 *
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *   tsv_filename will be the name of the TSV output
 *
 *********************************************************************/
void print_momo_results(int argc, char **argv, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  // Create directory for motifs to have a location
  create_directory(options);
  
  // Print HTML and TSV (must come first so it can set the E-values)
  FILE *html_file = fopen(options.html_path, "w");
  if (!html_file) {
    die("Couldn't open file %s for output.\n", options.html_path);
  }
  fclose(html_file);
  FILE *tsv_file = fopen(options.tsv_path, "w");
  if (!tsv_file) {
    die("Couldn't open file %s for output.\n", options.tsv_path);
  }
  print_momo_html_and_tsv_files(argc, argv, tsv_file, options, summary);

  // Print plain text.
  FILE *text_file = fopen(options.text_path, "w");
  if (!text_file) {
    die("Couldn't open file %s for output.\n", options.text_path);
  }
  print_momo_text_file(text_file, options, summary);
  fclose(text_file);
} // print_momo_results
