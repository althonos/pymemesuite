/********************************************************************
 * MOMO Portal
 ********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "display.h"
#include "dir.h"
#include "momo.h"
#include "momo-output.h"
#include "matrix.h"
#include "motif-in.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "binomial.h"
#include "meme.h"
#include "read_seq_file.h"
#include "momo-algorithm.h"
#include "momo-modl.h"

/* local variables */
#define DATA_HASH_SIZE 100003

/* print function for debugging */
void print_regexmotifs(ARRAYLST_T* regexmotifs,
                       int max_motifs,
                       MOMO_OPTIONS_T* options,
                       SUMMARY_T* summary) {
  int i;
  for (i = 0; i < arraylst_size(regexmotifs); ++i) {
    REGEX_MOTIF_T* regexmotif = arraylst_get(i,regexmotifs);
    print_matrix_to_terminal(regexmotif->residues, options, summary);
  }
}

/* print function for debugging */
void printOp(MODL_STEP_T step) {
  if (step.op == INIT_STEP) {
    printf("op == NULL_DL\n");
  } else if (step.op == MERGE) {
    printf("op == MERGE\n");
  } else if (step.op == ADD) {
    printf("op == ADD\n");
  } else if (step.op == APPEND) {
    printf("op == APPEND\n");
  } else if (step.op == REPLACE) {
    printf("op == REPLACE\n");
  } else {
    printf("op == REMOVE\n");
  }
}

/***********************************************************************
 Initialize regexmotif
 ***********************************************************************/
REGEX_MOTIF_T* init_regexmotif(MOMO_OPTIONS_T* options,
                               SUMMARY_T* summary) {
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  int i;
  REGEX_MOTIF_T* regexmotif = mm_malloc(sizeof(REGEX_MOTIF_T));
  regexmotif->conserved = (bool*) mm_malloc(sizeof(bool)*(width_no_center));
  for (i = 0; i < width_no_center; ++i) {
    regexmotif->conserved[i] = false;
  }
  regexmotif->residues = allocate_matrix(width_no_center, strlen(alph_letters));
  init_matrix(0.0, regexmotif->residues);
  
  return regexmotif;
}

/***********************************************************************
 Duplicate regexmotif
 ***********************************************************************/
REGEX_MOTIF_T* duplicate_regexmotif(REGEX_MOTIF_T* regexmotif,
                                    MOMO_OPTIONS_T* options,
                                    SUMMARY_T* summary) {
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  int i;
  REGEX_MOTIF_T* duplicate = mm_malloc(sizeof(REGEX_MOTIF_T));
  duplicate->conserved = (bool*) mm_malloc(sizeof(bool)*(width_no_center));
  for (i = 0; i < width_no_center; ++i) {
    duplicate->conserved[i] = regexmotif->conserved[i];
  }
  duplicate->residues = duplicate_matrix(regexmotif->residues);
  return duplicate;
}

/***********************************************************************
 Clean up regexmotif
 ***********************************************************************/
void cleanup_regexmotif(REGEX_MOTIF_T* regexmotif) {
  myfree(regexmotif->conserved);
  free_matrix(regexmotif->residues);
  myfree(regexmotif);
}

/***********************************************************************
 Initialize final regexmotif score
 ***********************************************************************/
FINAL_REGEX_MOTIF_SCORES_T* init_final_regexmotif_score(REGEX_MOTIF_T* motif, double score) {
  FINAL_REGEX_MOTIF_SCORES_T* final_regexmotif_score = mm_malloc(sizeof(FINAL_REGEX_MOTIF_SCORES_T));
  final_regexmotif_score->motif = motif;
  final_regexmotif_score->score = score;
  return final_regexmotif_score;
}

/***********************************************************************
 Clean up final regexmotif score
 ***********************************************************************/
void cleanup_final_regexmotif_score(FINAL_REGEX_MOTIF_SCORES_T* final_regexmotif_score) {
  myfree(final_regexmotif_score);
}

/***********************************************************************
 Initialize modl_step
 ***********************************************************************/
MODL_STEP_T* init_modl_step(MODL_OP_T op,
                            int mtfIdx1,
                            int mtfIdx2,
                            int position,
                            int residue,
                            double score) {
  MODL_STEP_T* modl_step = mm_malloc(sizeof(MODL_STEP_T));
  modl_step->op = op;
  modl_step->mtfIdx1 = mtfIdx1;
  modl_step->mtfIdx2 = mtfIdx2;
  modl_step->position = position;
  modl_step->residue = residue;
  modl_step->score = score;
  modl_step->n_tests = 0;
  return modl_step;
}

/***********************************************************************
 Clean up modl_step
 ***********************************************************************/
void cleanup_modl_step(MODL_STEP_T* modl_step) {
  myfree(modl_step);
}

/**
 * Returns the length of the residues for a given pos of a regexmotif.
 */
int get_residues_length(REGEX_MOTIF_T* regexmotif, SUMMARY_T* summary, int pos) {
  const char* alph_letters = summary->alph_letters;
  int i;
  int residues_length = 0;
  for (i = 0; i < strlen(alph_letters); ++i) {
    residues_length += (int) get_matrix_cell(pos, i, regexmotif->residues);
  }
  return residues_length;
}

/**
 * Computes the number of bits (description length) to encode the given regexmotif.
 */
double compute_regexmotif_numbits(REGEX_MOTIF_T* regexmotif,
                                  MOMO_OPTIONS_T* options,
                                  SUMMARY_T* summary) {
  int i;
  int j;
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  
  double result = 0.0;
  result += width_no_center; // number of bits to store conserved
  
  int num_conserved = 0; // number of bits to store residue representation
  for (i = 0; i < width_no_center; ++i) {
    if (regexmotif->conserved[i]) {
      num_conserved++;
    }
  }
  result += num_conserved;
  
  // figure out residue representation and number of bits to encode residues
  for (i = 0; i < width_no_center; ++i) {
    if (regexmotif->conserved[i]) {
      int residues_length = get_residues_length(regexmotif, summary, i);
      double binaryrep = (residues_length+1)*ceil(log(strlen(alph_letters))/log(2));
      double listrep = strlen(alph_letters);
      result += (binaryrep < listrep ? binaryrep : listrep);
    }
  }
  return result;
}

/**
 * Converts a regexmotif to a string representation.
 */
char* regexmotif_to_string(REGEX_MOTIF_T* regexmotif,
                           MOD_INFO_T* mod_info,
                           SUMMARY_T* summary,
                           MOMO_OPTIONS_T* options) {
  int i;
  int j;
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;

  // Required Length
  // the length of mod/center
  // 2 for _ surrounding mod/center
  // 1 for the '\0' character
  // 1 for non-conserved or single lettered positions
  // length + 2 for positions with multi residues (to encode [])
  int required_length = 3 + strlen(mod_info->mod_name);
  for (i=0; i < width_no_center; ++i) {
    int residues_length = (regexmotif->conserved[i]) ? get_residues_length(regexmotif, summary, i) : 0;
    if (residues_length == 0 || residues_length == 1) {
      required_length++;
    } else {
      required_length += residues_length + 2;
    }
  }
  
  // motif string creation
  char* result = mm_malloc(required_length);
  result[0] = '\0';
  for (i = 0; i < options->width; ++i) {
    int idx_no_center = (i < width_no_center / 2) ? i : i - 1;
    if (i == options->width / 2) {
      strncat(result, "_", 1);
      strncat(result, mod_info->mod_name, strlen(mod_info->mod_name));
      strncat(result, "_", 1);
    } else if (!regexmotif->conserved[idx_no_center]) {
      strncat(result, "x", 1);
    } else {
      int residues_length = get_residues_length(regexmotif, summary, idx_no_center);
      if (residues_length > 1) {
        strncat(result, "[", 1);
      }
      for (j = 0; j < strlen(alph_letters); ++j) {
        if (get_matrix_cell(idx_no_center, j, regexmotif->residues) == 1.0) {
          strncat(result, &alph_letters[j], 1);
        }
      }
      if (residues_length > 1) {
        strncat(result, "]", 1);
      }
    }
  }
  return result;
}

/**
 * Returns if the sequence contains the regexmotif. Returns false if regexmotif is null.
 */
bool sequence_contains_regexmotif(REGEX_MOTIF_T* regexmotif,
                                       char* sequence,
                                       MOMO_OPTIONS_T* options,
                                       SUMMARY_T* summary) {
  if (regexmotif == NULL) {
    return false;
  }
  int i;
  int width_no_center = options->width - 1;
  const char* alph_letters = summary->alph_letters;
  
  for (i = 0; i < width_no_center; ++i) {
    if (regexmotif->conserved[i]) {
      int char_idx = (i < width_no_center / 2) ? i : i+1;
      char* alph_letters_substring = strchr(alph_letters, sequence[char_idx]);
      
      if (alph_letters_substring) {
        int aa_idx = alph_letters_substring - alph_letters;
        if (get_matrix_cell(i, aa_idx, regexmotif->residues) != 1.0) {
          return false;
        }
      } else {
        return false;
      }
    }
  }
  return true;
}

/**
 * Computes the number of bits (description length) to encode all the peptides 
 * using the given list of motifs.
 */
double compute_peptides_dl(ARRAYLST_T* regexmotifs,
                           MATRIX_T* bg_freqs,
                           int max_motifs,
                           MOMO_OPTIONS_T* options,
                           SUMMARY_T* summary,
                           MOD_INFO_T* mod_info) {
  const char* alph_letters = summary->alph_letters;
  int i;
  int j;
  int k;
  int width_no_center = options->width - 1;
  
  // we want to compute a dl for each peptide
  double result = 0.0;
  
  // number of motifs
  double num_motifs = 0.0;
  for (i = 0; i < arraylst_size(regexmotifs); ++i) {
    num_motifs = num_motifs + 1.0;
  }
  
  // for each peptide
  for (i = 0; i < arraylst_size(mod_info->seq_list); ++i) {
    SEQ_T* seq_object = (SEQ_T*) ((options->eliminate_repeat_width) ? hash_get_entry_value(arraylst_get(i, mod_info->seq_list)) : arraylst_get(i, mod_info->seq_list));
    char* sequence = get_raw_sequence(seq_object);
    // number of bits to store does this sequence contain motifs?
    result += 1.0;
    
    // which of the motifs are in this sequence?
    bool has_motif = false;
    
    bool within_sequence[arraylst_size(regexmotifs)];
    for (j = 0; j < arraylst_size(regexmotifs); ++j) {
      within_sequence[j] = sequence_contains_regexmotif(arraylst_get(j, regexmotifs), sequence, options, summary);
      if (within_sequence[j]) {
        has_motif = true;
      }
    }
    
    result += has_motif ? num_motifs : 0;

    // cost of each residue
    for (j = 0; j < width_no_center; ++j) {
      int char_idx = (j < width_no_center / 2) ? j : j+1;
      int min_length_residues_list = 0;
      // get min_length_residues_list
      for (k = 0; k < arraylst_size(regexmotifs); ++k) {
        REGEX_MOTIF_T* regexmotif = arraylst_get(k,regexmotifs);
        if (within_sequence[k] && regexmotif->conserved[j]) {
          char* alph_letters_substring = strchr(alph_letters, sequence[char_idx]);
          if (alph_letters_substring) {
            int aa_idx = alph_letters_substring - alph_letters;
            
            if (get_matrix_cell(j, aa_idx, regexmotif->residues) == 1.0) {
              int residues_length = get_residues_length(regexmotif, summary, j);
              if (min_length_residues_list == 0 || residues_length < min_length_residues_list) {
                min_length_residues_list = residues_length;
              }
            }
          }
        }
      }
      // cost of residue using min_length_residue_list
      if (min_length_residues_list == 0) {
        if (strchr(alph_letters, sequence[char_idx])) {
          result -= log(get_matrix_cell(char_idx, strchr(alph_letters, sequence[char_idx]) - alph_letters, bg_freqs)) / log(2);
        }
      } else {
        result += log(min_length_residues_list) / log(2);
      }
    }
  }
  return result;
}

/**
 * Computes the number of bits (description length) to encode the given motifs and peptides.
 */
double compute_MDL(ARRAYLST_T* regexmotifs,
                   MATRIX_T* bg_freqs,
                   int max_motifs,
                   MOMO_OPTIONS_T* options,
                   SUMMARY_T* summary,
                   MOD_INFO_T* mod_info) {
  int i;
  double log2k = 7.0; // ceil log2(100) where k always = 100 for in MoDL implementation for some reason ?
  double motifs_dl = 0.0;
  for (i = 0; i < arraylst_size(regexmotifs); ++i) {
    motifs_dl += compute_regexmotif_numbits(arraylst_get(i, regexmotifs), options, summary);
  }
  double peptides_dl = compute_peptides_dl(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
  return log2k + motifs_dl + peptides_dl;
}

/**
 * Returns the resulting DL if the 'residue' was added at the 'position' 
 * to the specified motif at index 'number'. Modifies the list of motifs 
 * only if update is true.
 */
double add_residue(ARRAYLST_T* regexmotifs,
                   MATRIX_T* bg_freqs,
                   int max_motifs,
                   int number,
                   int position,
                   int residue,
                   MOMO_OPTIONS_T* options,
                   SUMMARY_T* summary,
                   MOD_INFO_T* mod_info,
                   bool update) {
  // Original state
  REGEX_MOTIF_T* regexmotif = arraylst_get(number, regexmotifs);
  double init_state = get_matrix_cell(position, residue, regexmotif->residues);
  bool init_conserved = regexmotif->conserved[position];
  
  regexmotif->conserved[position] = true;
  set_matrix_cell(position, residue, 1.0, regexmotif->residues);
  
  double dl = compute_MDL(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
  
  // Revert to original state if not update
  if (!update) {
    regexmotif->conserved[position] = init_conserved;
    set_matrix_cell(position, residue, init_state, regexmotif->residues);
  }
  return dl;
}

/**
 * Returns the resulting DL if the motif an index 'number' is removed. 
 * Modifies the list of motifs only if update is true.
 */
double remove_motif(ARRAYLST_T* regexmotifs,
                    MATRIX_T* bg_freqs,
                    int max_motifs,
                    int number,
                    MOMO_OPTIONS_T* options,
                    SUMMARY_T* summary,
                    MOD_INFO_T* mod_info,
                    bool update) {
  REGEX_MOTIF_T* regexmotif = arraylst_remove(number, regexmotifs);
  double dl = compute_MDL(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
  if (!update) {
    arraylst_put(number, regexmotif, regexmotifs);
  }
  return dl;
}

/**
 * Returns the resulting DL if the motifs an index 'mtfIdx1' and 'mtfIdx2' were removed 
 * and the merged result of the two motifs were appended to the end of the list. Modifies 
 * the list of motifs only if update is true.
 */
double merge_motifs(ARRAYLST_T* regexmotifs,
                    MATRIX_T* bg_freqs,
                    int max_motifs,
                    int mtfIdx1,
                    int mtfIdx2,
                    MOMO_OPTIONS_T* options,
                    SUMMARY_T* summary,
                    MOD_INFO_T* mod_info,
                    bool update) {

  int i;
  int j;
  
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  
  REGEX_MOTIF_T* regexmtf2 = arraylst_remove(mtfIdx2, regexmotifs);
  REGEX_MOTIF_T* regexmtf1 = arraylst_remove(mtfIdx1, regexmotifs);
  
  // generate a merged motif
  REGEX_MOTIF_T* mergedmotif = init_regexmotif(options, summary);
  for (i = 0; i < width_no_center; ++i) {
    mergedmotif->conserved[i] = (regexmtf1->conserved[i] || regexmtf2->conserved[i]);
    for (j = 0; j < strlen(alph_letters); ++j) {
      if (get_matrix_cell(i, j, regexmtf1->residues) == 1.0
          || get_matrix_cell(i, j, regexmtf2->residues) == 1.0) {
        set_matrix_cell(i, j, 1.0, mergedmotif->residues);
      }
    }
  }
  
  arraylst_add(mergedmotif, regexmotifs);
  
  double dl = compute_MDL(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
  
  if (update) {
    cleanup_regexmotif(regexmtf1);
    cleanup_regexmotif(regexmtf2);
  } else {
    arraylst_put(mtfIdx1, regexmtf1, regexmotifs);
    arraylst_put(mtfIdx2, regexmtf2, regexmotifs);
    arraylst_remove(arraylst_size(regexmotifs)-1, regexmotifs);
    cleanup_regexmotif(mergedmotif);
  }
  
  return dl;
}

/**
 * Returns the minimum DL (stored as a step) you can get by merging two motifs inside the list of motifs.
 */
MODL_STEP_T* get_best_merge(ARRAYLST_T* regexmotifs,
                            MATRIX_T* bg_freqs,
                            int max_motifs,
                            MOMO_OPTIONS_T* options,
                            SUMMARY_T* summary,
                            MOD_INFO_T* mod_info) {
  int i;
  int j;
  
  MODL_STEP_T* merge_step = init_modl_step(MERGE, 0, 0, 0, 0, INFINITY);
  for (i = 0; i < arraylst_size(regexmotifs); ++i) {
    for (j = i+1; j < arraylst_size(regexmotifs); ++j) {
      double dl = merge_motifs(regexmotifs, bg_freqs, max_motifs, i, j, options, summary, mod_info, false);
      merge_step->n_tests++;
      if (dl < merge_step->score) {
        merge_step->mtfIdx1 = i;
        merge_step->mtfIdx2 = j;
        merge_step->score = dl;
      }
    }
  }
  return merge_step;
}

/**
 * Returns the minimum DL (stored as a step) you can get by adding a single lettered motif to list of motifs.
 */
MODL_STEP_T* get_best_add(ARRAYLST_T* regexmotifs,
                          MATRIX_T* bg_freqs,
                          int max_motifs,
                          MOMO_OPTIONS_T* options,
                          SUMMARY_T* summary,
                          MOD_INFO_T* mod_info) {
  int i, j;
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  
  MODL_STEP_T* add_step = init_modl_step(ADD, 0, 0, 0, 0, INFINITY);
  
  if (arraylst_size(regexmotifs) < max_motifs) {
    REGEX_MOTIF_T* empty_motif = init_regexmotif(options, summary);
    arraylst_add(empty_motif, regexmotifs);
    
    for (i = 0; i < width_no_center; ++i) {
      for (j = 0; j < strlen(alph_letters); ++j) {
        double dl = add_residue(regexmotifs, bg_freqs, max_motifs,
                                arraylst_size(regexmotifs)-1, i, j, options, summary, mod_info, false);
	add_step->n_tests++;
        if (dl < add_step->score) {
          add_step->position = i;
          add_step->residue = j;
          add_step->score = dl;
        }
      }
    }
    
    arraylst_remove(arraylst_size(regexmotifs)-1, regexmotifs);
    cleanup_regexmotif(empty_motif);
  }
  
  return add_step;
}

/**
 * Returns the minimum DL (stored as a step) you can get by copying one of the motifs, 
 * then adding a single letter to the copy of the motif to the list of motifs. Also,
 * it allows the option of removing one of the motifs from the list other than the newly 
 * modified one. [See momo-modl.h struct modl_step] for more info.
 */
MODL_STEP_T* get_best_append(ARRAYLST_T* regexmotifs,
                             MATRIX_T* bg_freqs,
                             int max_motifs,
                             MOMO_OPTIONS_T* options,
                             SUMMARY_T* summary,
                             MOD_INFO_T* mod_info) {
  int i, j, k, l;
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  
  MODL_STEP_T* append_step = init_modl_step(APPEND, 0, 0, 0, 0, INFINITY);
  
  if (arraylst_size(regexmotifs) < max_motifs) {
    for (i = 0; i < arraylst_size(regexmotifs); ++i) {
      REGEX_MOTIF_T* dup = duplicate_regexmotif(arraylst_get(i, regexmotifs), options, summary);
      arraylst_add(dup, regexmotifs);
      for (j = 0; j < arraylst_size(regexmotifs) + 1; ++j) {
        REGEX_MOTIF_T* temp = NULL;
        if (j < arraylst_size(regexmotifs)-1) {
          temp = arraylst_remove(j, regexmotifs);
        }
        for (k = 0; k < width_no_center; ++k) {
          for (l = 0; l < strlen(alph_letters); ++l) {
            if (get_matrix_cell(k, l, dup->residues) != 1.0) {
              double dl = add_residue(regexmotifs, bg_freqs, max_motifs,
                                      arraylst_size(regexmotifs)-1, k, l, options, summary, mod_info, false);
	      append_step->n_tests++;
              if (dl < append_step->score) {
                append_step->mtfIdx1 = i;
                append_step->mtfIdx2 = j;
                append_step->position = k;
                append_step->residue = l;
                append_step->score = dl;
              }
            }
          }
        }
        if (temp) {
          arraylst_put(j, temp, regexmotifs);
        }
      }
      arraylst_remove(arraylst_size(regexmotifs)-1, regexmotifs);
      cleanup_regexmotif(dup);
    }
  }
  return append_step;
}

/**
 * Returns the minimum DL (stored as a step) you can get by adding a single letter 
 * to one of the motifs in the list of motifs (and moving it to the end). Also,
 * it allows the option of removing one of the motifs from the list other than the newly 
 * modified one. [See momo-modl.h struct modl_step] for more info.
 */
MODL_STEP_T* get_best_replace(ARRAYLST_T* regexmotifs,
                              MATRIX_T* bg_freqs,
                              int max_motifs,
                              MOMO_OPTIONS_T* options,
                              SUMMARY_T* summary,
                              MOD_INFO_T* mod_info) {
  int i;
  int j;
  int k;
  int l;
  const char* alph_letters = summary->alph_letters;
  int width_no_center = options->width - 1;
  
  MODL_STEP_T* replace_step = init_modl_step(REPLACE, 0, 0, 0, 0, INFINITY);
  
  for (i = 0; i < arraylst_size(regexmotifs); ++i) {
    REGEX_MOTIF_T* replace_motif = arraylst_remove(i, regexmotifs);
    arraylst_add(replace_motif, regexmotifs);
    for (j = 0; j < arraylst_size(regexmotifs); ++j) {
      REGEX_MOTIF_T* temp = NULL;
      if (j < arraylst_size(regexmotifs)-1) {
        temp = arraylst_remove(j, regexmotifs);
      }
      for (k = 0; k < width_no_center; ++k) {
        for (l = 0; l < strlen(alph_letters); ++l) {
          if (get_matrix_cell(k, l, replace_motif->residues) != 1.0) {
            double dl = add_residue(regexmotifs, bg_freqs, max_motifs,
                                    arraylst_size(regexmotifs)-1, k, l, options, summary, mod_info, false);
	    replace_step->n_tests++;
            if (dl < replace_step->score) {
              replace_step->mtfIdx1 = i;
              replace_step->mtfIdx2 = j;
              replace_step->position = k;
              replace_step->residue = l;
              replace_step->score = dl;
            }
          }
        }
      }
      if (temp) {
        arraylst_put(j, temp, regexmotifs);
      }
    }
    arraylst_remove(arraylst_size(regexmotifs)-1, regexmotifs);
    arraylst_put(i, replace_motif, regexmotifs);
  }

  return replace_step;
}

/**
 * Returns the minimum DL (stored as a step) you can get by removing one of the motifs in the list of motifs.
 */
MODL_STEP_T* get_best_remove(ARRAYLST_T* regexmotifs,
                             MATRIX_T* bg_freqs,
                             int max_motifs,
                             MOMO_OPTIONS_T* options,
                             SUMMARY_T* summary,
                             MOD_INFO_T* mod_info) {
  int i;
  MODL_STEP_T* remove_step = init_modl_step(REMOVE, 0, 0, 0, 0, INFINITY);
  for (i = 0; i < arraylst_size(regexmotifs) - 1; ++i) {
    double dl = remove_motif(regexmotifs, bg_freqs, max_motifs, i, options, summary, mod_info, false);
    remove_step->n_tests++;
    if (dl < remove_step->score) {
      remove_step->mtfIdx1 = i;
      remove_step->score = dl;
    }
  }
  return remove_step;
}

/**
 * Modifies list of motifs by executing the given step.
 */
void do_step(MODL_STEP_T* step,
             ARRAYLST_T* regexmotifs,
             MATRIX_T* bg_freqs,
             int max_motifs,
             MOMO_OPTIONS_T* options,
             SUMMARY_T* summary,
             MOD_INFO_T* mod_info) {
  
  if (step->op == MERGE) {
    merge_motifs(regexmotifs, bg_freqs, max_motifs, step->mtfIdx1, step->mtfIdx2,
                 options, summary, mod_info, true);
  } else if (step->op == ADD) {
    REGEX_MOTIF_T* empty_motif = init_regexmotif(options, summary);
    arraylst_add(empty_motif, regexmotifs);
    add_residue(regexmotifs, bg_freqs, max_motifs, arraylst_size(regexmotifs)-1,
                step->position, step->residue, options, summary, mod_info, true);
  } else if (step->op == APPEND) {
    REGEX_MOTIF_T* dup = duplicate_regexmotif(arraylst_get(step->mtfIdx1, regexmotifs), options, summary);
    arraylst_add(dup, regexmotifs);
    add_residue(regexmotifs, bg_freqs, max_motifs, arraylst_size(regexmotifs)-1,
                step->position, step->residue, options, summary, mod_info, true);
    if (step->mtfIdx2 < arraylst_size(regexmotifs)-1) {
      arraylst_remove(step->mtfIdx2, regexmotifs);
    }
  } else if (step->op == REPLACE) {
    REGEX_MOTIF_T* replace_motif = arraylst_remove(step->mtfIdx1, regexmotifs);
    arraylst_add(replace_motif, regexmotifs);
    add_residue(regexmotifs, bg_freqs, max_motifs, arraylst_size(regexmotifs)-1,
                step->position, step->residue, options, summary, mod_info, true);
    if (step->mtfIdx2 < arraylst_size(regexmotifs)-1) {
      arraylst_remove(step->mtfIdx2, regexmotifs);
    }
  } else if (step->op == REMOVE) {
    arraylst_remove(step->mtfIdx1, regexmotifs);
  }
}

/**
 * Figure out the best set of motifs and convert this into motifinfos.
 * Motifs are ordered by the increase in DL if they are removed.
 */
ARRAYLST_T* get_best_motif_list(ARRAYLST_T* best_ops,
                                double minDL,
                                MATRIX_T* bg_freqs,
                                int max_motifs,
                                SUMMARY_T* summary,
                                MOMO_OPTIONS_T* options,
                                MOD_INFO_T * mod_info) {
  int i;
  int j;
  int k;
  
  ARRAYLST_T* seq_list = mod_info->seq_list;
  
  ARRAYLST_T* finalmotifs = arraylst_create();
  ARRAYLST_T* temp = arraylst_create();
  int most_affected = 0;
  for (i = 0; i < arraylst_size(best_ops); ++i) {
    // Apply the step
    MODL_STEP_T* step = arraylst_get(i, best_ops);
    do_step(step, temp, bg_freqs, max_motifs, options, summary, mod_info);
    
    // Update final motifs only if score = minDL and has more affected sequences
    if (step->score == minDL) {
      // Get number of sequences affected by this set
      int num_affected = 0;
      for (j = 0; j < arraylst_size(seq_list); ++j) {
        SEQ_T* seq_object = (SEQ_T*) ((options->eliminate_repeat_width) ? hash_get_entry_value(arraylst_get(j, seq_list)) : arraylst_get(j, seq_list));
        char* sequence = get_raw_sequence(seq_object);
        bool contains_motif = false;
        for (k = 0; k < arraylst_size(temp); ++k) {
          if (sequence_contains_regexmotif(arraylst_get(k, temp), sequence, options, summary)) {
            contains_motif = true;
          }
        }
        if (contains_motif) {
          num_affected++;
        }
      }
      // Update if number affected > previously
      if (num_affected > most_affected) {
        // clear the sequences in final motifs
        for (j = 0; j < arraylst_size(finalmotifs); ++j) {
          REGEX_MOTIF_T* motif = arraylst_get(j, finalmotifs);
          cleanup_regexmotif(motif);
        }
        arraylst_clear(NULL, finalmotifs);
        // duplicate temp into final motifs
        for (j = 0; j < arraylst_size(temp); ++j) {
          REGEX_MOTIF_T* motif = arraylst_get(j, temp);
          arraylst_add(duplicate_regexmotif(motif, options, summary), finalmotifs);
        }
      }
    }
  }
  // Cleanup temp
  for (i = 0; i < arraylst_size(temp); ++i) {
    REGEX_MOTIF_T* motif = arraylst_get(i, temp);
    cleanup_regexmotif(motif);
  }
  arraylst_destroy(NULL, temp);
  
  return finalmotifs;
}

/*
 * final regexmotif comparator
 * This compares 2 pointers to FINAL_REGEX_MOTIF_SCORES_T*.
 */
static int final_regex_motif_score_compare(const void *mtf1_p, const void *mtf2_p) {
  FINAL_REGEX_MOTIF_SCORES_T *mtf1, *mtf2;
  mtf1 = *((FINAL_REGEX_MOTIF_SCORES_T**) mtf1_p);
  mtf2 = *((FINAL_REGEX_MOTIF_SCORES_T**) mtf2_p);
  // ascending: 1 - 2, descending: 2 - 1
  return (mtf2->score - mtf1->score);
}

/**
 * Figure out the best set of motifs and convert this into motifinfos.
 * Motifs are ordered by the increase in DL if they are removed.
 */
void set_motifinfo_list(ARRAYLST_T* finalmotifs,
                        double minDL,
                        MATRIX_T* bg_freqs,
                        int max_motifs,
                        SUMMARY_T* summary,
                        MOMO_OPTIONS_T* options,
                        MOD_INFO_T * mod_info,
			int n_tests) {
  int i;
  int j;
  
  // We would like to rank each final motif by the increase in DL when they are removed.
  // Larger increases are better.
  ARRAYLST_T* finalscores = arraylst_create();;
  for (i = 0; i < arraylst_size(finalmotifs); ++i) {
    // temporarily remove this motif
    REGEX_MOTIF_T* regexmotif = arraylst_remove(i, finalmotifs);
    // compute the score
    double removeDL = compute_MDL(finalmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    double score = removeDL - minDL;
    FINAL_REGEX_MOTIF_SCORES_T* final_regexmotif_score = init_final_regexmotif_score(regexmotif, score);
    arraylst_add(final_regexmotif_score, finalscores);
    // put the motif back
    arraylst_put(i, regexmotif, finalmotifs);
  }
  
  // sort the motifs
  arraylst_qsort(final_regex_motif_score_compare, finalscores);
  
  ARRAYLST_T* seq_list = mod_info->seq_list;
  ARRAYLST_T* bg_seq_list = mod_info->bg_seq_list;
  
  // Each seq should only be added to one motif. True if seq is part of a motif.
  bool status[arraylst_size(seq_list)];
  bool bg_status[arraylst_size(bg_seq_list)];
  for (i = 0; i < arraylst_size(seq_list); ++i) status[i] = false;
  for (i = 0; i < arraylst_size(bg_seq_list); ++i) bg_status[i] = false;
  int fg_matches_total=0, bg_matches_total=0;
  
  // Loop over motifs.
  for (i = 0; i < arraylst_size(finalscores); ++i) {
    FINAL_REGEX_MOTIF_SCORES_T* final_regexmotif_score = arraylst_get(i, finalscores);
    REGEX_MOTIF_T* motif = final_regexmotif_score->motif;
    double score = final_regexmotif_score->score;
    int bg_matches, fg_matches;
    
    char* mod_name = regexmotif_to_string(motif, mod_info, summary, options);
    ARRAYLST_T* matching_seqs = arraylst_create();
    ARRAYLST_T* matching_chars = arraylst_create();
    fg_matches = 0;
    for (j = 0; j < arraylst_size(seq_list); ++j) {
      SEQ_T* seq_object = (SEQ_T*) ((options->eliminate_repeat_width) ? hash_get_entry_value(arraylst_get(j, seq_list)) : arraylst_get(j, seq_list));
      char* sequence = get_raw_sequence(seq_object);
      if (!status[j] && sequence_contains_regexmotif(motif, sequence, options, summary)) {
        arraylst_add(arraylst_get(j, seq_list), matching_seqs);
        arraylst_add(sequence, matching_chars);
        fg_matches++;
        fg_matches_total++;		// matches of SET of motifs 
        status[j] = true;
      }
    }
    bg_matches = 0;
    for (j = 0; j < arraylst_size(bg_seq_list); ++j) {
      SEQ_T* bg_seq_object = (SEQ_T*) ((options->eliminate_repeat_width) ? hash_get_entry_value(arraylst_get(j, bg_seq_list)) : arraylst_get(j, bg_seq_list));
      char* bg_sequence = get_raw_sequence(bg_seq_object);
      if (!bg_status[j] && sequence_contains_regexmotif(motif, bg_sequence, options, summary)) {
        bg_matches++;
        bg_matches_total++;		// matches of SET of motifs 
        bg_status[j] = true;
      }
    }
    
    // make the motifinfo
    const char* alph_letters = summary->alph_letters;
    
    // Create the frequency matrix
    MATRIX_T* freqs = NULL;
    freqs = get_count_matrix(freqs, matching_seqs, NULL, options, summary);
    normalize_rows(0.0, freqs);
    
    arraylst_destroy(NULL, matching_seqs);
    
    // Create the motif
    MOTIF_INFO_T* motifinfo = mm_malloc(sizeof(MOTIF_INFO_T));
    motifinfo->motif = allocate_motif(mod_name, "", summary->alph, freqs, NULL);
    motifinfo->seqs = matching_chars;
    motifinfo->fg_matches = arraylst_size(matching_chars);
    motifinfo->fg_size = arraylst_size(mod_info->seq_list);
    motifinfo->bg_matches = bg_matches;
    motifinfo->bg_size = arraylst_size(mod_info->bg_seq_list);
    motifinfo->score = score;
    motifinfo->n_tests = n_tests;
    // Get the counts of all matches in foreground sequences.
    motifinfo->afg_size = arraylst_size(mod_info->seq_list);
    motifinfo->afg_matches = get_counts_from_motifid(mod_name, mod_info->seq_list, options);
    // Get the counts of all matches in background sequences.
    motifinfo->abg_size = arraylst_size(mod_info->bg_seq_list);
    motifinfo->abg_matches = get_counts_from_motifid(mod_name, mod_info->bg_seq_list, options);

    arraylst_add(motifinfo, mod_info->motifinfos);

    // clean up
    free_matrix(freqs);
  }

  // Add total matches to the SET of motifs to motifinfos.
  //for (i = 0; i < arraylst_size(mod_info->motifinfos); ++i) {
  //  MOTIF_INFO_T* motifinfo = arraylst_get(i, mod_info->motifinfos);
  //  motifinfo->afg_matches = fg_matches_total;
  //  motifinfo->afg_size = motifinfo->fg_size;
  //  motifinfo->abg_matches = bg_matches_total;
  //  motifinfo->abg_size = motifinfo->bg_size;
  //}
  
  // Clean up
  for (i = 0; i < arraylst_size(finalscores); ++i) {
    FINAL_REGEX_MOTIF_SCORES_T* final_regexmotif_score = arraylst_get(i, finalscores);
    cleanup_final_regexmotif_score(final_regexmotif_score);
  }
  arraylst_destroy(NULL, finalscores);
}

/**
 * Creates and stores the motifs found using the modl
 * algorithm for a given mod.
 */
void create_modl_motifs(SUMMARY_T* summary,
                        MOMO_OPTIONS_T* options,
                        MOD_INFO_T * mod_info
                        ) {

  const char* alph_letters = summary->alph_letters;
  int i, j, k;
  int max_motifs = options->max_motifs;
  int max_iterations = options->max_iterations;
  int max_no_decrease = options->max_no_decrease;
  
  // background frequencies
  MATRIX_T* bg_freqs = NULL;
  bg_freqs = get_count_matrix(bg_freqs, mod_info->bg_seq_list, NULL, options, summary);
  int alphsize = strlen(alph_letters);
  double prior = 0.05;
  for (i = 0; i < options->width; ++i ) {
    for (j = 0; j < alphsize; ++j) {
      // add prior
      set_matrix_cell(i,j,(get_matrix_cell(i,j,bg_freqs)+prior)/(arraylst_size(mod_info->bg_seq_list)+alphsize*prior),bg_freqs);
    }
  }
  
  int no_decrease = 0;
  
  ARRAYLST_T* regexmotifs = arraylst_create();
  ARRAYLST_T* best_ops = arraylst_create();
  
  double minDL = compute_MDL(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
  MODL_STEP_T* init_step = init_modl_step(INIT_STEP, 0, 0, 0, 0, minDL);
  arraylst_add(init_step, best_ops);
  
  // NOTE: The MoDL implementation from the Ritz Paper does something funny (potential bug?)
  // no decrease is checked after every time a op is added to best_ops, but the check
  // for size of ops < max iterations is only checked once everytime in the loop...
  // NOTE2: Each round through the loop has 3 potential places to add an op... This is
  // inconsistent with what is described by the paper, which seems to imply you can only have
  // one op per round through the loop.... It also makes this while loop
  // unnecessarily long and complicated...
  //   1. merging (optional)
  //   2. add, append(+remove), replace(+remove)
  //   3. remove (optional)
  // NOTE3: The possible operations are different than as described in the paper...
  // NOTE4: I dont understand why score of last element == minDL even a condition for loop terminiation ??
  int n_tests = 0;
  while (arraylst_size(best_ops)-1 <= max_iterations &&
         (((MODL_STEP_T*) arraylst_peek(best_ops))->score == minDL || no_decrease < max_no_decrease)) {
    // Try merge existing motifs
    MODL_STEP_T* merge_step = get_best_merge(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    n_tests += merge_step->n_tests;
    // Do the merge if it decreases the DL (compared to the last DL)
    MODL_STEP_T* last_op = arraylst_peek(best_ops);
    if (merge_step->score < last_op->score) {
      arraylst_add(merge_step, best_ops);
      do_step(merge_step, regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
      if (merge_step->score > minDL) { // yes... I do realize that >= makes more sense... but apparently this is what modl does ?
        no_decrease++;
      } else {
        no_decrease = 0;
        minDL = merge_step->score;
      }
    }
    
    // Add, append, replace [append and replace have options to remove as well]
    MODL_STEP_T* add_step = get_best_add(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    n_tests += add_step->n_tests;
    MODL_STEP_T* append_step = get_best_append(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    n_tests += append_step->n_tests;
    MODL_STEP_T* replace_step = get_best_replace(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    n_tests += replace_step->n_tests;

    // Figure out which is the best of the three and apply.
    MODL_STEP_T* best_add_append_replace = add_step;
    if (arraylst_size(regexmotifs) > 0) {
      if (append_step->score < best_add_append_replace->score) {
        best_add_append_replace = append_step;
      }
      if (replace_step->score < best_add_append_replace->score) {
        best_add_append_replace = replace_step;
      }
    }
    
    // Perform the best_add_append_replace regardless of whether or not there is a decrease in DL...
    if (best_add_append_replace->score > minDL) {
      no_decrease++;
    } else {
      no_decrease = 0;
      minDL = best_add_append_replace->score;
    }
    arraylst_add(best_add_append_replace, best_ops);
    do_step(best_add_append_replace, regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    
    // Try remove a motif if it helps decrease the DL...
    MODL_STEP_T* remove_step = get_best_remove(regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
    n_tests += remove_step->n_tests;
    // Do the remove if it decreases the DL (compared to the last DL)
    last_op = arraylst_peek(best_ops);
    if (remove_step->score < best_add_append_replace->score) {
      arraylst_add(remove_step, best_ops);
      do_step(remove_step, regexmotifs, bg_freqs, max_motifs, options, summary, mod_info);
      if (remove_step->score > minDL) {
        no_decrease++;
      } else {
        no_decrease = 0;
        minDL = remove_step->score;
      }
    }
  }
  
  mod_info->modl_ops = best_ops;
  
  ARRAYLST_T* finalmotifs = get_best_motif_list(best_ops, minDL, bg_freqs, max_motifs, summary, options, mod_info);
  set_motifinfo_list(finalmotifs, minDL, bg_freqs, max_motifs, summary, options, mod_info, n_tests);
  
  // cleanup
  for (i = 0; i < arraylst_size(regexmotifs); ++i) {
    REGEX_MOTIF_T* regexmotif = arraylst_get(i, regexmotifs);
    cleanup_regexmotif(regexmotif);
  }
  arraylst_destroy(NULL, regexmotifs);
  for (i = 0; i < arraylst_size(finalmotifs); ++i) {
    REGEX_MOTIF_T* regexmotif = arraylst_get(i, finalmotifs);
    cleanup_regexmotif(regexmotif);
  }
  arraylst_destroy(NULL, finalmotifs);
  free_matrix(bg_freqs);

}
