/********************************************************************
 * FILE: gomo.c
 * AUTHOR: Fabian Buske
 * CREATE DATE: 18/06/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, UQ
 *
 * GOMo is an implementation of the algorithm described in
 * "Associating transcription factor binding site motifs with target
 * Go terms and target genes"
 * authors: Mikael Boden and Timothy L. Bailey
 * published: Nucl. Acids Res (2008)
 *
 ********************************************************************/

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "mtwist.h"
// matrix.h has some strange requirement that it be defined before arrays.h
// as one of the .h files below includes arrays.h (I have no idea which) 
// matrix.h needs to stay at the top...
#include "matrix.h"
#include "alphabet.h"
#include "array-list.h"
#include "binary-search.h"
#include "buffer.h"
#include "ceqlogo.h"
#include "cisml.h"
#include "ceqlogo.h"
#include "config.h"
#include "dir.h"
#include "fasta-io.h"
#include "gomo_highlight.h"
#include "hash_table.h"
#include "io.h"
#include "linked-list.h"
#include "merger.h"
#include "motif-in.h"
#include "projrel.h"
#include "qvalue.h"
#include "ranksum_test.h"
#include "read_csv.h"
#include "red-black-tree.h"
#include "simple-getopt.h"
#include "string-builder.h"
#include "string-list.h"
#include "xml-util.h"
#include "xml-out.h"

#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a<b)?b:a
#define bool2txt(b)   (b) ?"true":"false"

/*********************************************************************************
 * Constants 
 ********************************************************************************/
#define GO_BLOCK_SIZE 100
#define DEFAULT_Q_VALUE_THRESHOLD 0.05  /* Q-value threshold */
#define BIN_ALLOC_CHUNK 10 /* number of tuples to increment bins by */
#define DEFAULT_SHUFFLES 1000
#define DEFAULT_MIN_GENE_COUNT 1
const char *program_name = "gomo"; /* the program name */
static char *default_out_dir = "gomo_out";  /* default name of output */
static const char *RES_DIRNAME = "res";
static const char *XML_FILENAME = "gomo.xml";
static const char *TSV_FILENAME = "gomo.tsv";
static const char *HTML_STYLESHEET = "gomo-to-html.xsl";
static const char *HTML_FILENAME = "gomo.html";
static const char *LOGO_PREFIX = "logo";
static const char *LOGO_SUFFIX = ".png";
//for get_gomap_path, which is not used. 
/*********************************************************************************
 * Globals
 ********************************************************************************/
VERBOSE_T verbosity = NORMAL_VERBOSE;
bool status = true;
time_t status_last = 0; //last time a status message was output
time_t status_delay = 5; //minimum time between status messages

#define DISPLAY_STATUS(status_msg_format, ...) { \
  if (status) { \
    fprintf(stderr, status_msg_format, __VA_ARGS__); \
    status_last = time(NULL); \
  } \
} \

#define SKIPABLE_STATUS(status_msg_format, ...) { \
  if (status) { \
    time_t status_new = time(NULL); \
    if (status_new > (status_last + status_delay)) { \
      fprintf(stderr, status_msg_format, __VA_ARGS__); \
      status_last = status_new; \
    } \
  } \
} 

/*************************************************************************
 * gomo_dts()
 *************************************************************************/
static char *gomo_dts = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
    "<!-- Document definition -->\n"
    "<!DOCTYPE gomo[\n"
    "\t<!ELEMENT gomo (program, motif*)>\n"
    "\t<!ATTLIST gomo xmlns:xsi CDATA #IMPLIED\n"
    "\t\t\t version CDATA #REQUIRED\n"
    "\t\t\t release CDATA #REQUIRED\n"
    "\t\t\t>\n"
    "\t\t<!ELEMENT program (gomapfile, seqscorefile*)>\n"
    "\t\t<!ATTLIST program\n"
    "\t\t\t\t name CDATA #REQUIRED\n"
    "\t\t\t\t cmd CDATA #REQUIRED\n"
    "\t\t\t\t gene_url CDATA #IMPLIED\n"
    "\t\t\t\t outdir CDATA #REQUIRED\n"
    "\t\t\t\t clobber CDATA #REQUIRED\n"
    "\t\t\t\t text_only CDATA #REQUIRED\n"
    "\t\t\t\t use_e_values CDATA #REQUIRED\n"
    "\t\t\t\t score_e_thresh CDATA #REQUIRED\n"
    "\t\t\t\t min_gene_count CDATA #REQUIRED\n"
    "\t\t\t\t motifs CDATA #IMPLIED\n"
    "\t\t\t\t shuffle_scores CDATA #REQUIRED\n"
    "\t\t\t\t q_threshold CDATA #REQUIRED\n"
    "\t\t\t\t>\n"
    "\t\t\t<!ELEMENT gomapfile EMPTY>\n"
    "\t\t\t<!ATTLIST gomapfile path CDATA #REQUIRED >\n"
    "\t\t\t<!ELEMENT seqscorefile EMPTY>\n"
    "\t\t\t<!ATTLIST seqscorefile path CDATA #REQUIRED >\n"
    "\t\t<!ELEMENT motif (goterm*)>\n"
    "\t\t<!ATTLIST motif\n"
    "\t\t\t\t id CDATA #REQUIRED\n"
    "\t\t\t\t genecount CDATA #REQUIRED\n"
    "\t\t\t\t logo CDATA #IMPLIED\n"
    "\t\t\t\t>\n"
    "\t\t\t<!ELEMENT goterm (gene*)>\n"
    "\t\t\t<!ATTLIST goterm\n"
    "\t\t\t\t\t id CDATA #REQUIRED\n"
    "\t\t\t\t\t score CDATA #REQUIRED\n"
    "\t\t\t\t\t pvalue CDATA #REQUIRED\n"
    "\t\t\t\t\t qvalue CDATA #REQUIRED\n"
    "\t\t\t\t\t annotated CDATA #REQUIRED\n"
    "\t\t\t\t\t group CDATA #REQUIRED\n"
    "\t\t\t\t\t nabove CDATA #REQUIRED\n"
    "\t\t\t\t\t nbelow CDATA #REQUIRED\n"
    "\t\t\t\t\t implied (u|y|n) #REQUIRED\n"
    "\t\t\t\t\t name CDATA #REQUIRED\n"
    "\t\t\t\t\t>\n"
    "\t\t\t\t<!ELEMENT gene EMPTY>\n"
    "\t\t\t\t<!ATTLIST gene\n"
    "\t\t\t\t\t\t id CDATA #REQUIRED\n"
    "\t\t\t\t\t\t rank CDATA #REQUIRED\n"
    "\t\t\t\t\t\t>\n"
    "]>\n";

/***********************************************************************************
 * Data types
 **********************************************************************************/
/*
 * seq mapping
 * Provides a mapping from one sequence to another
 * allowing shuffling of labels on sequence scores.
 */
typedef struct seq_mapping SHUFFLER_T;
struct seq_mapping {
  int *map;           // list that contains a mapping for all loaded sequences
                      // provides a conversion from one sequence to another
                      // a sequential numbering(map[i] = i) provides the
                      // identity mapping.
  int map_length;     // the number of entries in map
  int *rank;          // list that contains lookup for go map sequences
                      // generated using the seq_map to map a gomap sequence
                      // to a rank in the sorted score list
  int rank_length;    // the number of entries in rank
};

/*
 * ama score
 * Store the score of a motif for a sequence. This score is loaded from a cisml file
 * and converted to an e-value if it was a p-value, or the gene score is used.
 * At the time of loading a unique integer identifier is generated for each
 * sequence identifier (aka gene/operon ID) in a one to one mapping, this is 
 * known as the load number of the sequence. 
 */
typedef struct ama_score AMA_SCORE_T;
struct ama_score {
  char *seq_id; // the sequence identifier
  int seq_load_num; // the load number of the sequence 
  double score; // the score for a sequence, it could be an evalue or a score 
};

/*
 * ama list
 * Stores a list of ama scores of sequences for a single motif. The scores 
 * are sorted best to worst
 */
typedef struct ama_list AMA_LIST_T;
struct ama_list {
  AMA_SCORE_T *list;
  int length;
};

typedef struct ama_group AMA_GROUP_T;
struct ama_group {
  MOTIF_T *motif;
  AMA_LIST_T **species;
  int species_count;
};

/* 
 * gomo result
 * Used for collating the result just before it is output.
 */
typedef struct gomo_result GOMO_RESULT_T;
struct gomo_result {
  int term_load_num; // index in go map
  char* term_id; // go term 
  double gomo_score; // geometric mean of rank sum p-values
  double empirical_pvalue; // score for go term
};

/*
 * sequence list
 * A list of sequences which have been annotated with a GO term.
 * The list uses the sequence load numbers.
 */
typedef struct seq_list_t SEQ_LIST_T;
struct seq_list_t {
  int length; // the number of annotated sequences
  int *seq_load_nums; //the load numbers of the annotated sequences
};

/*
 * A tuple contains the multi-species GOMo score associated with 
 * a given GO term and the counts of Null model items that have been 
 * better. Note: the go_id is not copied on creation so it
 * should not be deallocated on destruction.
 */
typedef struct tuple_t TUPLE_T;
struct tuple_t {
  int term_load_num; //the index in the go map
  char *term_id;                // the GO term ID
  double gomo_score;    // the gomo score of the term 
  int counts; // number of counts
};

/*
 * Each GO term has a number of genes associated with it, called
 * its "gene_count".  The BIN_T is designed to contain
 * all of the tuples (go index, goid, gomo score, counts) entries for
 * all GO terms within the given gene count range.  The tuples are
 * stored in an array that always has an allocated size that is a 
 * multiple of BIN_ALLOC_CHUNK. 
 */
typedef struct bin_t BIN_T;
struct bin_t {
  int low_gene_count; // low bound on gene count (inclusive)
  int high_gene_count; // high bound on gene count (inclusive)
  TUPLE_T *tuples;     // array of tuples within gene count range 
  int tuple_count;  // the number of valid tuples in tuples
  int total_counts; // the total number of counts added to the bin
};

/*
 * The list of bins sorted in order of range ascending. The bins will 
 * have non-overlapping ranges and there will be a bin for any possible
 * gene count. 
 */
typedef struct bins_t BINS_T;
struct bins_t {
  BIN_T **values; 
  int length; //allocated length of values array
  int total_tuples; //total number of tuples in bins
};


/************************************************************************
 * Bins, Bin and Tuple structure modification methods
 ***********************************************************************/

/*
 * Comparison function for two structures of TUPLE_T.
 * Compares based on the gomo score of the tuple.
 */
int tuple_pvalue_compar(const void *v1, const void *v2) 
{
  TUPLE_T *t1 = (TUPLE_T*)v1;
  TUPLE_T *t2 = (TUPLE_T*)v2;
  if (t1->gomo_score < t2->gomo_score) {
    return -1;
  } else if (t1->gomo_score > t2->gomo_score) {
    return 1;
  } else {
    return 0;
  }
}

/*
 * Creates a BIN_T for holding tuples for GO terms that have been annotated
 * with a number of genes between the inclusive range of 'low_gene_count' and
 * 'high_gene_count'.
 *
 * 1. Allocate memory for the bin.
 * 2. Set all the values of the bin.
 * 3. Return the bin.
 */
BIN_T* bin_create(int low_gene_count, int high_gene_count) 
{
  assert(low_gene_count < high_gene_count);
  BIN_T *bin = mm_malloc(sizeof(BIN_T));
  bin->low_gene_count = low_gene_count;
  bin->high_gene_count = high_gene_count;
  bin->tuple_count = 0;
  bin->tuples = NULL;
  bin->total_counts = 0;
  return bin;
}

/*
 * Destroys a BIN_T created by bin_create.
 *
 * 1. If there are any tuples, free them.
 * 2. Free the bin.
 */
void bin_destroy(BIN_T *bin) 
{
  assert(bin != NULL);
  if (bin->tuple_count > 0) free(bin->tuples);
  free(bin);
}

/*
 * Adds a tuple created from the values 'go_index', 'go_term' and 'score'
 * to the specified bin with a counts value of zero.
 *
 * 1. Ensure there is enough space in the bin->tuples to add a tuple.
 * 2. Set the values of the tuple.
 * 3. Increment the tuple count.
 */
void bin_add_tuple(int term_load_num, char *term_id, double score, BIN_T* bin) 
{
  if (bin->tuple_count % BIN_ALLOC_CHUNK == 0) {
    int new_size = bin->tuple_count + BIN_ALLOC_CHUNK;
    Resize(bin->tuples, new_size, TUPLE_T);
  }
  TUPLE_T* tuple = (bin->tuples)+(bin->tuple_count);
  tuple->term_load_num = term_load_num;
  tuple->term_id = term_id;
  tuple->gomo_score = score;
  tuple->counts = 0;
  bin->tuple_count += 1;
}

/*
 * Compares a key with an entry. The key is expected to be a pointer
 * to a double representing a gomo score. The entry is expected to be a
 * pointer to a TUPLE_T which also contains a gomo score. The two gomo scores
 * are compared and the function returns -1, 0 or 1 respectively if the
 * key's gomo score is smaller, equal to, or greater, than the entry's gomo 
 * score.
 */
int bin_bsearch_compar(const void *key, const void *entry) {
  double score = *((double*)key);
  TUPLE_T *tuple = (TUPLE_T*)entry;
  if (score < tuple->gomo_score) {
    return -1;
  } else if (score == tuple->gomo_score) {
    return 0;
  }
  return 1;
}

/*
 * Adds a count for a score generated by the null model. The counts allow 
 * emperical calculation of the probability of a score occuring by random 
 * chance (p-value). This assumes the bin has been sorted.
 *
 * 1. binary search the tuples to find a tuple that has the same score or 
 *    a score just larger.
 * 2. if such a tuple exists then increment its count.
 */
void bin_add_count(double score, BIN_T* bin)
{
  int pos = binary_search(&score, bin->tuples, bin->tuple_count, 
      sizeof(TUPLE_T), bin_bsearch_compar);
  if (pos < 0) pos = -(pos + 1);
  if (pos < bin->tuple_count) {
    (bin->tuples+pos)->counts += 1;
  }
  bin->total_counts += 1;
}


/*
 * Sort the tuples in the bin by p-value ascending.
 */
void bin_sort(BIN_T *bin) 
{
  qsort(bin->tuples, bin->tuple_count, sizeof(TUPLE_T), tuple_pvalue_compar);
}

/*
 * Create a bins data structure.
 * 
 * 1. Allocate memory for the bins structure and for the pointers to the 
 *    bin types.
 * 2. Create the bin(s) so that there is a bin starting at zero and at 
 *    each division so that the entire possible range of gene counts has 
 *    a bin allocated.
 * 3. Return the bins structure.
 */
BINS_T* bins_create(int divisions_length, // the number of places to split 
    int *divisions                         // places to spit the bins
    ) 
{
  assert(divisions_length >= 0);
  assert(divisions_length == 0 || divisions != NULL);
  BINS_T *bins = mm_malloc(sizeof(BINS_T));
  bins->length = 1 + divisions_length; 
  bins->values = mm_malloc(sizeof(BIN_T*)*bins->length);
  //create bins upto but not including the divisions
  int i, low;
  for (low = 0, i = 0; i < divisions_length; ++i) {
    assert(low < divisions[i]);
    bins->values[i] = bin_create(low, divisions[i]-1);
    low = divisions[i];
  }
  //create final bin
  bins->values[i] = bin_create(low, INT_MAX);
  bins->total_tuples = 0;
  return bins;
}

/*
 * Destroy a bin structure created by bins_create.
 */
void bins_destroy(BINS_T *bins) 
{
  assert(bins != NULL);
  int i;
  for (i = 0; i < bins->length; i++) {
    bin_destroy(bins->values[i]);
  }
  free(bins->values);
  free(bins);
}

/*
 * Compares a gene count (key) with a bin (value) to determine
 * if the gene count would fit in the bin.
 */
int bins_bsearch_compar(const void *key, const void *entry) {
  int gene_count = *((int*)key);
  BIN_T *bin = *(BIN_T**)entry;
  if (gene_count >= bin->low_gene_count) {
    if (gene_count <= bin->high_gene_count) {
      return 0;
    } else {
      return 1;
    }
  }
  return -1;
}

/*
 * Use to get a bin for a gene_count. If the bin
 * doesn't exist then it returns NULL.
 */
BIN_T* bins_get_bin(int gene_count, BINS_T *bins) 
{
  assert(bins != NULL);
  assert(gene_count >= 0);
  int index = binary_search(&gene_count, bins->values, bins->length, sizeof(BIN_T*), bins_bsearch_compar);
  if (index < 0) return NULL;
  return bins->values[index];
}

/*
 * Bins a tuple based on the number of genes that are associated with the 'go_term'.
 */
void bins_add_bin_tuple(int go_index, char *go_term, double score, int gene_count, BINS_T *bins) 
{
  assert(go_term != NULL);
  assert(bins != NULL);
  assert(gene_count >= 0); 
  bin_add_tuple(go_index, go_term, score, bins_get_bin(gene_count, bins));
  bins->total_tuples += 1;
}

/*
 * Bins a count based on the number of genes that are associated with the 'go_term'.
 * This assumes that the bins have been sorted.
 */
void bins_add_bin_count(double score, int gene_count, BINS_T *bins) 
{
  bin_add_count(score, bins_get_bin(gene_count, bins));
}

/*
 * Sorts all the bin structures held by bins.
 * This must be called before counts can be added to bins.
 */
void bins_sort(BINS_T *bins) 
{
  assert(bins != NULL);
  int i;
  for (i = 0; i < bins->length; ++i) {
    bin_sort(bins->values[i]);
  }
}

/*****************************************************************************
 * Methods relating to the SHUFFLER_T 
 *
 ****************************************************************************/
/*
 * shuffler_create
 * Creates a shuffler. The shuffling needs to be defined before this is used.
 */
SHUFFLER_T *shuffler_create(int gomap_seq_count, int total_seq_count) 
{
  assert(total_seq_count >= gomap_seq_count);
  SHUFFLER_T *shuffler = (SHUFFLER_T*)mm_malloc(sizeof(SHUFFLER_T));
  shuffler->map = (int*)mm_malloc(sizeof(int)*total_seq_count);
  shuffler->map_length = total_seq_count;
  shuffler->rank = (int*)mm_malloc(sizeof(int)*gomap_seq_count);
  shuffler->rank_length = gomap_seq_count;
  return shuffler;
}

/*
 * shuffler_destroy
 * Destroys a shuffler.
 */
void shuffler_destroy(SHUFFLER_T *shuffler) 
{
  free(shuffler->map);
  free(shuffler->rank);
}

/*
 * shuffler_identity
 * Sets the shuffler to the identity mapping
 */
void shuffler_identity(SHUFFLER_T *shuffler) 
{
  int i;
  for (i = 0; i < shuffler->map_length; ++i) {
    shuffler->map[i] = i;
  }
}

/*
 * shuffler_shuffle
 * Shuffles the map. Assumes that all map entries
 * are unique from a previous call to shuffler_identity
 */
void shuffler_shuffle(SHUFFLER_T *shuffler) 
{
  int i, s, *list, length;
  list = shuffler->map;
  length = shuffler->map_length;
  for (i = 0; i < length; ++i) {
    s = (int)(mt_ldrand() * (length - i)) + i;
    int tmp = list[i];
    list[i] = list[s];
    list[s] = tmp;
  }
}

/*************************************************************************
 * gomo_result_cmp()
 *
 * Compares two gomo results with respect to its scores
 *************************************************************************/
int gomo_result_cmp(const void *v1, const void *v2) 
{
  assert(v1 != NULL);
  assert(v2 != NULL);
  GOMO_RESULT_T s1 = *(GOMO_RESULT_T *)v1;
  GOMO_RESULT_T s2 = *(GOMO_RESULT_T *)v2;
  if (s1.empirical_pvalue > s2.empirical_pvalue){
    return 1;
  } else if (s1.empirical_pvalue < s2.empirical_pvalue){
    return -1;
  } else {
    return 0;
  }
}

/*************************************************************************
 * start_xml_output()
 *
 * initiate the output; returns command line
 *************************************************************************/
char *start_xml_output(FILE *xml_output,
    uint32_t seed,
    int argc,
    char *argv[],
    char *go_map_file,
    ARRAYLST_T *ama_paths,
    char *outdir,
    bool clobber,
    bool text_only,
    bool use_e_values,
    double score_e_thresh,
    int min_gene_count,
    ARRAYLST_T *motifs,
    int shuffle_scores,
    double q_threshold,
    char *gene_url
    )
{
  fputs(gomo_dts, xml_output);
  fputs("<gomo xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
      "version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n", xml_output);

  STR_T *buf = str_create(10);
  // create command line
  char *command = NULL;
  int lc = strlen(program_name);
  int alloc_size = lc + 1;
  int i;
  for (i=1;i<argc;++i){
    alloc_size += strlen(argv[i])+1;
  }
  command = mm_malloc(sizeof(char)*alloc_size);
  strncpy(command, program_name, lc);
  char* ptr;
  if (argc > 0) {
    for (i = 1; i < argc; ++i) {
      command[lc++] = ' ';
      ptr = argv[i];
      while (*ptr != '\0') {
        command[lc++] = *(ptr++);
      }
    }
  }
  command[lc] = '\0';
  // output program state 
  fprintf(xml_output,"\t<program name=\"%s\"\n\t\t\t seed=\"%" PRIu32 "\"\n",
      xmlify(program_name, buf, true), seed);
  fprintf(xml_output, "\t\t\t cmd=\"%s\"\n", xmlify(command, buf, true));
  fprintf(xml_output, "\t\t\t gene_url=\"%s\"\n", xmlify(gene_url, buf, true));
  fprintf(xml_output,
      "\t\t\t outdir=\"%s\" clobber=\"%s\" text_only=\"%s\"\n"
      "\t\t\t use_e_values=\"%s\" score_e_thresh=\"%g\"\n"
      "\t\t\t min_gene_count=\"%d\"", 
      xmlify(outdir, buf, true), bool2txt(clobber), bool2txt(text_only), 
      bool2txt(use_e_values), score_e_thresh, min_gene_count);
  if (motifs != NULL) {
    fprintf(xml_output, " motifs=\"%s", xmlify((char*)arraylst_get(0, motifs), buf, true));
    for (i = 1; i < arraylst_size(motifs); ++i) {
      fprintf(xml_output, " %s", xmlify((char*)arraylst_get(i, motifs), buf, true));
    }
    fprintf(xml_output, "\"");
  }
  fprintf(xml_output, "\n\t\t\t shuffle_scores=\"%d\" q_threshold=\"%g\"\n\t\t\t>\n", shuffle_scores, q_threshold);
  fprintf(xml_output,"\t\t<gomapfile path=\"%s\" />\n", xmlify(go_map_file, buf, true));
  for (i = 0; i < arraylst_size(ama_paths); ++i) {
    fprintf(xml_output,"\t\t<seqscorefile path=\"%s\" />\n", xmlify((char*)arraylst_get(i, ama_paths), buf, true));
  }
  fprintf(xml_output, "\t</program>\n");
  // clean up
  //free(command);
  str_destroy(buf, false);
  return(command);
}

double get_ama_score(const void *scores_ptr, int index) {
  AMA_LIST_T *scores = (AMA_LIST_T*)scores_ptr;
  return scores->list[index].score; 
}
bool get_false(const void *ignore1, int ignore2) {
  return false; 
}

/*
 * ama_evalue_compar
 * Compare two AMA evalues. Used for sorting smallest to biggest
 * which is best to worst. If the first is smaller, the same,
 * or larger than it will return repectively -1, 0, 1.
 */
int ama_evalue_compar(const void *v1, const void *v2) {
  AMA_SCORE_T *s1 = (AMA_SCORE_T*)v1;
  AMA_SCORE_T *s2 = (AMA_SCORE_T*)v2;
  if (s1->score < s2->score) {
    return -1;
  } else if (s1->score > s2->score) {
    return 1;
  }
  return 0;
}

/*
 * ama_evalue_compar
 * Compare two AMA scores. Used for sorting biggest to smallest
 * which is best to worst. If the first is larger, the same,
 * or smaller than it will return repectively -1, 0, 1.
 */
int ama_score_compar(const void *v1, const void *v2) {
  AMA_SCORE_T *s1 = (AMA_SCORE_T*)v1;
  AMA_SCORE_T *s2 = (AMA_SCORE_T*)v2;
  if (s1->score > s2->score) {
    return -1;
  } else if (s1->score < s2->score) {
    return 1;
  }
  return 0;
}

/*
 * create_ama_score_list
 * Makes a copy of only the sequence scores and the identifying
 * sequence load number (the order in which the sequence was loaded from file).
 * Additionally sorts the sequence list by the score best to worst
 */
AMA_LIST_T* create_ama_score_list(
    PATTERN_T *pattern,     // pattern to copy data from 
    bool use_e_values, // choose between the score source
    double ama_e_threshold, // if the score source is converted p-values
                            // then this is the threshold at which the e-values
                            // are set to the maximum (meaning no association)
    RBTREE_T *seq_ids       // the already loaded sequence ids
    ) 
{
  AMA_LIST_T *ama_scores;
  AMA_SCORE_T *input;
  SCANNED_SEQUENCE_T *sseq;
  int i, seq_count, seq_loading_num;
  bool created;
  RBNODE_T *node;
  
  //find out how many scored sequences there are 
  seq_count = get_pattern_num_scanned_sequences(pattern);
  //set the sequence loading number to the next one
  seq_loading_num = rbtree_size(seq_ids);
  //allocate memory for the ama scores 
  ama_scores = mm_malloc(sizeof(AMA_LIST_T));
  ama_scores->list = mm_malloc(sizeof(AMA_SCORE_T)*seq_count);
  ama_scores->length = seq_count;
  //loop over all sequence scores
  for (i = 0, input = ama_scores->list; i < seq_count; ++i, ++input) {
    sseq = get_pattern_scanned_sequences(pattern)[i];
    //make a copy of the sequence name (as I made the decision to not provide a copy function to the constructor)
    char * sseq_name;
    copy_string(&sseq_name, get_scanned_sequence_name(sseq));
    //lookup the seq id
    node = rbtree_lookup(seq_ids, sseq_name, true, &created);
    if (created) {// assign it a loading number
      rbtree_set(seq_ids, node, &seq_loading_num);
      ++seq_loading_num;
    } else {
      free(sseq_name);
    }
    //set the load number for this ama score
    input->seq_id = (char*)rbtree_key(node);
    input->seq_load_num = *((int*)rbtree_value(node));
    // select between copying e-values (calculated from pvalues) and sequence scores
    if (use_e_values) {
      if (!has_scanned_sequence_pvalue(sseq)) 
        die("The sequence score file does not provide a "
            "p-value for sequence %s (pattern %s). "
            "Use --gs to prompt GOMo working on "
            "gene scores rather than gene-score p-values.\n",
            get_scanned_sequence_name(sseq),
            get_pattern_accession(pattern));
      // combine the p-value with the sequence count to make an e-value
      input->score = get_scanned_sequence_pvalue(sseq)*seq_count;
      // if the e-value is above the specified threshold then consider the
      // association with the motif to be by random chance and score it accordingly
      if (ama_e_threshold > 0 && input->score > ama_e_threshold) {
        input->score = seq_count;
      }
    } else {
      if (!has_scanned_sequence_score(sseq)) 
        die("The sequence score file does not provide a "
            "score for sequence %s (pattern %s). "
            "Remove --gs to prompt GOMo working on "
            "gene-score p-values rather than gene scores.\n",
            get_scanned_sequence_name(sseq),
            get_pattern_accession(pattern));
      input->score = get_scanned_sequence_score(sseq);
    }
  }
  // sort the inputs in order of best to worst
  if (use_e_values) {
    //for evalues smallest is best
    qsort(ama_scores->list, ama_scores->length, sizeof(AMA_SCORE_T), ama_evalue_compar);
  } else {
    //for scores biggest is best
    qsort(ama_scores->list, ama_scores->length, sizeof(AMA_SCORE_T), ama_score_compar);
  }
  return ama_scores;
}

/*
 * destroy_ama_score_list
 * Frees all memory associated with the passed AMA_LIST_T.
 */
void destroy_ama_score_list(AMA_LIST_T *ama_scores) {
  free(ama_scores->list);
  free(ama_scores);
}

/*
 * create_ama_score_group
 * Allocates memory for an ama score group
 * the score lists are left linking to nothing
 */
AMA_GROUP_T* create_ama_score_group(int species_count) {
  AMA_GROUP_T *group;

  group = (AMA_GROUP_T*)mm_malloc(sizeof(AMA_GROUP_T));
  group->motif = NULL;
  group->species_count = species_count;
  group->species = (AMA_LIST_T**)mm_calloc(species_count, sizeof(AMA_LIST_T*));

  return group;
}

/*
 * destroy_ama_score_group
 * Frees all memory associated with the passed AMA_GROUP_T.
 * It takes a void pointer to be compatible with the red black
 * tree.
 */
void destroy_ama_score_group(void *v1) {
  int i;
  AMA_GROUP_T *group;
  
  group = (AMA_GROUP_T*)v1;
  for (i = 0; i < group->species_count; ++i) {
    if (group->species[i] == NULL) break;
    destroy_ama_score_list(group->species[i]);
  }
  free(group->species);
  if (group->motif) destroy_motif(group->motif);
  free(group);
}


/*****************************************************************************
 * File reading helper functions
 ****************************************************************************/

inline static int is_newline(void *config, int letter) {
  return (letter == '\n' || letter == '\r');
}

inline static int is_tab_or_nl(void *config, int letter) {
  return (letter == '\t' || letter == '\n' || letter == '\r');
}

/*
 * read a specified number of bytes
 */
char* read_chunk(BUF_T *buffer, FILE *fp, int len) {
  int i, c;
  char *chunk;

  if (len <= 0) return NULL;
  chunk = mm_malloc(sizeof(char)*(len+1));
  for (i = 0; i < len; ++i) {
    while ((c = buf_getc(buffer)) == -1) {
      if (feof(fp)) {
        die("File ends with %d bytes still expected\n", (len - i));
      }
      buf_compact(buffer);
      buf_fread(buffer, fp);
      if (ferror(fp)) {
          die("Error occured while reading the file, error was given as: %s\n", strerror(ferror(fp)));
      }
      buf_flip(buffer);
    }
    chunk[i] = (char)c; //narrowing conversion
  }
  chunk[i] = '\0';
  return chunk;
}

/*
 * skip seperators
 * keep skipping until encounting something that isn't a newline or a tab
 * returns true if encountering a new line, false if we only find tabs
 */
bool skip_seperators(BUF_T *buffer, FILE *fp) {
  bool newline;
  int c;
  newline = false;
  while (true) {
    while ((c = buf_getc(buffer)) != -1) {
      if (c == '\n' || c == '\r') {
        newline = true;
      } else if (c != '\t') {
        //put back non seperator
        buf_set_position(buffer, buf_position(buffer) - 1);
        return newline;
      }
    }
    if (feof(fp)) {
      return newline;
    }
    buf_compact(buffer);
    buf_fread(buffer, fp);
    if (ferror(fp)) {
      die("Error while reading, error given as %s\n", strerror(ferror(fp)));
    }
    buf_flip(buffer);
  }
}

/******************************************************************************
 * load_map()
 *
 * Reads in the go map for a gomo run
 * Reads the url for looking up information about gene ids.
 * Reads in unique ids
 *****************************************************************************/
SEQ_LIST_T* load_map(
    char* go_map,       // map to load IN
    int min_gene_count, // minimum gene count to be considered for loading IN
    char** gene_url,    // url at which the gene ids can be looked up OUT
    RBTREE_T **go_ids,  // go ids mapped to their load number OUT
    RBTREE_T **seq_ids  // seq ids mapped to their load number OUT
    ) 
{
  assert(go_map != NULL);
  assert(min_gene_count > 0);

  BUF_T *buffer;
  FILE *fp;
  int i, term_loading_num, seq_loading_num, size, *seq_load_num;
  char *token;
  ARRAYLST_T *term_seq_load_nums;
  bool created, newline;
  SEQ_LIST_T *map;
  RBNODE_T *term_node, *seq_node;

  *go_ids = rbtree_create(rbtree_strcmp,NULL,free,rbtree_intcpy,free);
  *seq_ids = rbtree_create(rbtree_strcasecmp,NULL,free,rbtree_intcpy,free);
  term_seq_load_nums = arraylst_create();

  fp = fopen(go_map, "r");
  if (fp == NULL) {
    die("Failed to open file \"%s\", error was given as %s\n", go_map, strerror(errno));
  }
  buffer = buf_create(100);
  buf_flip(buffer);

  DISPLAY_STATUS("Loading map from \"%s\".\n", go_map);

  map = NULL;
  term_node = NULL;
  term_loading_num = 0;
  seq_loading_num = 0;
  *gene_url = buf_fread_token(buffer, fp, is_newline, NULL, false, NULL, 0, NULL);
  if (*gene_url == NULL) {
    die("Error while reading file \"%s\", error was given as %s\n", go_map, strerror(ferror(fp)));
  }
  newline = skip_seperators(buffer,fp);
  while (!feof(fp) || buf_remaining(buffer)) {
    token = buf_fread_token(buffer, fp, is_tab_or_nl, NULL, false, NULL, 0, NULL);
    if (token == NULL) {
      die("Error while reading file \"%s\", error was given as %s\n", go_map, strerror(ferror(fp)));
    }
    if (newline) {//new line was previous char
      // process previous term
      if (term_loading_num > 0) {
        if (arraylst_size(term_seq_load_nums) < min_gene_count) { //this term has too few annotations so we don't want it
          rbtree_delete(*go_ids, term_node, NULL, NULL);
          --term_loading_num;
        } else {
          size = (term_loading_num-1);
          if (size % GO_BLOCK_SIZE == 0) {
            Resize(map, size + GO_BLOCK_SIZE, SEQ_LIST_T);
          }
          map[size].length = arraylst_size(term_seq_load_nums);
          map[size].seq_load_nums = mm_malloc(sizeof(int)*arraylst_size(term_seq_load_nums));
          for (i = 0, seq_load_num = map[size].seq_load_nums; i < arraylst_size(term_seq_load_nums); ++i, ++seq_load_num) {
            *seq_load_num = *((int*)arraylst_get(i, term_seq_load_nums));
          }
        }
        arraylst_clear(NULL, term_seq_load_nums);
      }
      //lookup the term id
      term_node = rbtree_lookup(*go_ids, token, true, &created);
      if (!created) {
        die("Term '%s' appears more than once!", token);
      }
      //assign a loading number
      rbtree_set(*go_ids, term_node, &term_loading_num);
      ++term_loading_num;
    } else { //tab was previous char
      //lookup the seq id
      seq_node = rbtree_lookup(*seq_ids, token, true, &created);
      if (created) {// assign it a loading number
        rbtree_set(*seq_ids, seq_node, &seq_loading_num);
        ++seq_loading_num;
      } else { // destroy the unused token
        free(token);
      }
      //keep track of the sequences for this go term
      arraylst_add(rbtree_value(seq_node), term_seq_load_nums);
    }
    newline = skip_seperators(buffer,fp);
  }
  // process final term
  if (term_loading_num > 0) {
    size = (term_loading_num-1);
    Resize(map, term_loading_num, SEQ_LIST_T);
    map[size].length = arraylst_size(term_seq_load_nums);
    map[size].seq_load_nums = mm_malloc(sizeof(int)*arraylst_size(term_seq_load_nums));
    for (i = 0, seq_load_num = map[size].seq_load_nums; i < arraylst_size(term_seq_load_nums); ++i, ++seq_load_num) {
      *seq_load_num = *((int*)arraylst_get(i, term_seq_load_nums));
    }
  }

  //cleanup
  fclose(fp);
  arraylst_destroy(NULL, term_seq_load_nums);

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Loaded %d unique GO terms (IDs) with %d or more annotations and %d annotating sequence IDs.\n", 
        rbtree_size(*go_ids), min_gene_count, rbtree_size(*seq_ids));
  }

  return map;
}


/******************************************************************************
 * load_ama_scores
 * Loads CisML files containing average motif afinity scores for sequences and
 * extracts the sequence scores, the type of which is dependent on the use_e_values
 * parameter. Creates a hash table with the motif as key and the sequence scores 
 * for each species as value. If filter_motifs is not null then it will only store
 * those motifs.
 *****************************************************************************/
void load_ama_scores(
    ARRAYLST_T *ama_paths,      // the list of score files to source ama scores
    bool use_e_values,     // the type of score to extract from the CisML file
    double ama_e_threshold,     // the threshold for considering random for a e-value
    ARRAYLST_T *filter_motifs,  // the motifs to limit analysis to, or null
    RBTREE_T *seq_ids,          // the already loaded sequence ids
    RBTREE_T ** motifs          // the data loaded from the motifs
  ) 
{
  assert(ama_paths != NULL);
  assert(arraylst_size(ama_paths) > 0);
  assert(seq_ids != NULL);
  assert(motifs != NULL);
  // some useful variables
  int i, j;
  PATTERN_T *pattern;
  char *accession, *ama_file;
  AMA_GROUP_T *m_scores;
  RBNODE_T *node;
  bool created;
  CISML_T *file;
  // load shared motifs 
  *motifs = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, destroy_ama_score_group);
  for (i = 0; i < arraylst_size(ama_paths); ++i) {
    ama_file = (char*)arraylst_get(i, ama_paths);
    DISPLAY_STATUS("Loading scores from \"%s\".\n", ama_file);
    file = read_cisml(ama_file);
    if (file == NULL) {
      // Wasn't cisml XML
      die("The scoring file needs to be in CisML format: %s",
          ama_file);
    }
    if (verbosity >= HIGHER_VERBOSE) {
      fprintf(stderr, "CisML file '%s' for species %d contains %d motifs with sequence scores.\n", 
          ama_file, i, get_cisml_num_patterns(file));
    }
    if (i == 0) {
      for (j = 0; j < get_cisml_num_patterns(file); ++j) {
        pattern = get_cisml_patterns(file)[j];  
        accession = get_pattern_accession(pattern); // motif identifier
        if (filter_motifs != NULL) {
          // if this motif hasn't been selected than skip it
          if (arraylst_bsearch(accession, arraylst_compar_txt, filter_motifs) < 0) continue;
        }
        // allocate memory to for the per-species lists of ama scores for this motif
        m_scores = create_ama_score_group(arraylst_size(ama_paths));
        // make a copy of everything useful in the pattern
        m_scores->species[0] = create_ama_score_list(pattern, use_e_values, ama_e_threshold, seq_ids);
        // associate the per-species sequence scores with the motif
        node = rbtree_lookup(*motifs, accession, true, &created);
        if (!created) {
          // the motif already has scores, this is strange so we should warn the user
          fprintf(stderr, "Pattern %s appears more than once in CisML file '%s'."
              " Only the first occurance will be tested.\n", accession,
              ama_file);
          // clean up
          destroy_ama_score_group(m_scores);
        } else {
          rbtree_set(*motifs, node, m_scores);
        }
      }
    } else {
      for (j = 0; j < get_cisml_num_patterns(file); ++j) {
        pattern = get_cisml_patterns(file)[j];  
        accession = get_pattern_accession(pattern);
        // check the motif has been in previous cisml files, otherwise skip
        if ((node = rbtree_lookup(*motifs, accession, false, NULL)) != NULL) {
          // we've seen this motif before, which is good
          m_scores = rbtree_value(node);
          if (m_scores->species[i-1] == NULL) {
            // the previous cisml file has not seen this motif though so we must remove it
            // as we can only work on motifs that all species have data for 
            rbtree_delete(*motifs, node, NULL, NULL);
            continue;
          }
          if (m_scores->species[i] != NULL) {
            // the motif already has scores, this is strange so we should warn the user
            fprintf(stderr, "Pattern %s appears more than once in CisML file '%s'."
                " Only the first occurance will be tested.\n", accession,
                ama_file);
          } else {
            // make a copy of everything useful in the pattern
            m_scores->species[i] = create_ama_score_list(pattern, use_e_values, ama_e_threshold, seq_ids);
          }
        }
      }
    }
    // clean up
    free_cisml(file);
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Loaded %d motifs which have sequence scores for all species.\n", 
        rbtree_size(*motifs));
  }
}

/*************************************************************************
 * calculate_ranksum_pvalues
 * Calculates a p-value for each GO-term using the Mann-Whitney U test.
 *
 * For each GO-term the sequences are classified into two groupings, those
 * that are annotated with the GO-term and those that are not. The ranksum
 * test (aka Mann-Whitney U test) is run to determine if the scores from the 
 * two classifications are from the same distribution giving a p-value.
 *
 *************************************************************************/
double *calculate_ranksum_pvalues(
    RBTREE_T *go_ids,
    SEQ_LIST_T *go_map,
    AMA_LIST_T *seq_scores,
    SHUFFLER_T *shuffler,
    bool null_score
  ) 
{
  double ta_obs;
  int i, j, class_a, seq_count;
  RBNODE_T *node;
  
  seq_count = seq_scores->length;

  // create mapping of sequence num to its index in the ordered list of scores
  for (i = 0; i < shuffler->rank_length; ++i) {
    shuffler->rank[i] = -1;
  }
  for (i = 0; i < seq_scores->length; ++i) {
    int index = seq_scores->list[i].seq_load_num;
    index = shuffler->map[index];
    if (index < shuffler->rank_length) {
      shuffler->rank[index] = i;
    }
  }

  // Copy the scores to the special object for ranksum routines.  
  // Calculate the (adjusted) rank of each score.
  // By adjusted rank this means that all consecutive equal scores
  // will have their new ranks made the average of their ranks.
  // Note that the group is not used, but it is set to false anyway.
  RSD_T *ranksum_dataset = get_ranksum_dataset2(seq_scores, get_ama_score, NULL, get_false, seq_count);

  int skipped = 0;
  double *go_scores = (double*)mm_malloc(sizeof(double)*rbtree_size(go_ids));
  // for each go term mark all the mapped sequences as being in class a 
  for (node = rbtree_first(go_ids); node != NULL; node = rbtree_next(node)) {
    //get the term loading number
    i = *((int*)rbtree_value(node));
    //initilize the stats we need to do a rank-sum from stats 
    class_a = 0; // number of elements in class a
    ta_obs = 0.0;       // sum of ranks in class a

    // look up the sequence ids associated with the go id 
    SEQ_LIST_T *seqids = go_map+i;

    // for each sequence id look up its rank
    for (j = 0; j < seqids->length; ++j) {
      // lookup the position of the sequence id and
      // extract the adjusted rank to sum for the rank-sum test
      // also keep count of how many we've summed
      int index = shuffler->rank[seqids->seq_load_nums[j]];
      if (index != -1) {
        ta_obs += get_ranksum_rank(ranksum_dataset, index);
        ++class_a;
      }
    }
    // cannot calculate if all are in one classification 
    if (class_a < seq_count && class_a > 0) {
      RSR_T *ranksum_result = ranksum_from_stats(seq_count, class_a, ta_obs);
      go_scores[i] = exp(RSR_get_log_p_left(ranksum_result));
      destroy_rsr(ranksum_result);
    } else {
      //no data
      if (null_score) {
        // as we're doing lots of samples a random function should average out
        go_scores[i] = mt_ldrand();
      } else {
        // only a few samples so can't use rand.
        go_scores[i] = 0.5; 
      }
      ++skipped;
      if (verbosity >= HIGHER_VERBOSE) {
        char *go_term = (char*)rbtree_key(node);
        fprintf(stderr, "GO term %s skipped as either all sequences were in its class or none were (%d / %d).\n", go_term, class_a, seq_count);
      }
    }
  }
  if (verbosity >= HIGH_VERBOSE) {
    if (skipped > 0)
      fprintf(stderr,"%d GO terms have assigned genes which lack a sequence entry"
          " or for other reasons have no information to work with. These are "
          "skipped!\n",skipped);
  }
  
  // clean up
  destroy_rsd(ranksum_dataset);

  return go_scores;
}

/*************************************************************************
 * calculate_scores_and_bin
 * Compute the GOMo multiple-species score for the motif for each GO-term.
 * Those scores are saved in "bins" indexed by number of genes annotated
 * to a GO-term.
 *
 *   1) Loops over species calling calculate_ranksum_pvalues, which 
 *      returns an array with all the single-species GOMo scores.
 *   2) Computes the multiple-species GOMo score for each GOterm by
 *          computing the geometric mean over species.
 *   3) If the data is not from the null model then tuples are added
 *      recording the goid and score, otherwise counts are added to the
 *      tuples 
 ************************************************************************/
void calculate_scores_and_bin(
    RBTREE_T *go_ids, 
    SEQ_LIST_T* go_map, 
    AMA_GROUP_T *m_scores, //all the scores for the motif
    bool is_nmod,
    SHUFFLER_T *shuffler,
    BINS_T *bins 
  ) 
{
  int go_load_num, species_index;
  char *go_id;
  RBNODE_T *node;
  // allocate memory for storing a p value for every go term per species
  double ** go_scores =  (double**)mm_malloc(sizeof(double*)*(m_scores->species_count));
  // loop on species
  for (species_index = 0; species_index < m_scores->species_count; ++species_index) {
    AMA_LIST_T *seq_scores = m_scores->species[species_index];
    go_scores[species_index] = seq_scores ? 
      calculate_ranksum_pvalues(go_ids, go_map, seq_scores, shuffler, is_nmod) : NULL;
  }
  // loop on go ids
  for (node = rbtree_first(go_ids); node != NULL; node = rbtree_next(node)) {
    //get the term and the load number
    go_id = (char*)rbtree_key(node);
    go_load_num = *((int*)rbtree_value(node));
    //calculate geometric mean
    double score = 1.0;
    int count = 0;
    for (species_index = 0; species_index < m_scores->species_count && score > 0; ++species_index) {
      //Note: if their was no information to calculate the p-value it will have been estimated to be 
      //0.5 if it was a real sample or 
      //drand_mt if it was a null model sample (which will average out to 0.5 over many tests)
      if (go_scores[species_index]) {
        score *= go_scores[species_index][go_load_num];
        count++;
      }
    }
    //score = pow(score, 1.0 / m_scores->species_count);// take nth root
    score = count>1 ? pow(score, 1.0 / count) : score;// take nth root
    //add tuple
    int gene_count = go_map[go_load_num].length;
    if (is_nmod) {
      bins_add_bin_count(score, gene_count, bins);
    } else {
      bins_add_bin_tuple(go_load_num, go_id, score, gene_count, bins);
    }
  }
  // clean up
  for (species_index = 0; species_index < m_scores->species_count; ++species_index) {
    if (go_scores[species_index]) free(go_scores[species_index]);
  }
  free(go_scores);
}

/*
 * integer comparison function
 */
int intcmp(const void *v1, const void *v2) {
  int i1, i2;
  i1 = *((int*)v1);
  i2 = *((int*)v2);
  if (i1 < i2) {
    return -1;
  } else if (i1 == i2) {
    return 0;
  } else {
    return 1;
  }
}

/**************************************************************
 * calculate_qvalues_and_output
 * Current motif p-values are computed empirically for each gene-count bin,
 * and outputs all GOterm GOMo scores (etc) for the motif. 
 * 
 * The empirical p-values for all GO-terms in a bin are estimated
 * using the counts of null model p-values that scored lower or equal. 
 * Then the p-values of all bins are input to the compute_q_value() function.
 * The E-values are computed by multiplying the empirical p-value by the
 * number of p-values.
 *************************************************************/
void calculate_qvalues_and_output(
    char *motif_id, 
    char *logo_path,
    BINS_T *bins, 
    FILE *txt_output,
    FILE *xml_output,
    bool text_only,
    double q_threshold,
    SEQ_LIST_T* go_map, // map
    AMA_LIST_T *seq_scores, // sequence scores for the last species
    int best_genes_max // the maximum number of best genes to select, may have less selected
    ) 
{
  assert(motif_id != NULL);
  assert(bins != NULL);
  
  STR_T *buf = str_create(10);
  GOMO_RESULT_T *results = (GOMO_RESULT_T*)mm_malloc(sizeof(GOMO_RESULT_T)*bins->total_tuples);

  int i, j, k, num_dupes, running_total; 
  GOMO_RESULT_T *result = results;
  for (i = 0; i < bins->length; ++i) {
    BIN_T *bin = bins->values[i];
    running_total = 0;
    for (j = 0; j < bin->tuple_count; j += num_dupes) {
      num_dupes = 0;
      do {
        running_total += bin->tuples[j+num_dupes].counts;
        ++num_dupes;
      } while (j+num_dupes < bin->tuple_count && bin->tuples[j].gomo_score == bin->tuples[j+num_dupes].gomo_score);
      double empirical_pvalue  = (double)(max(running_total, 1)) / (double) (bin->total_counts + 1);
      for (k = 0; k < num_dupes; ++k) { 
        result->term_load_num = bin->tuples[j+k].term_load_num;
        result->term_id = bin->tuples[j+k].term_id;
        result->gomo_score = bin->tuples[j].gomo_score;
        result->empirical_pvalue = empirical_pvalue; 
        ++result;
      }
    }
  }
  qsort(results, bins->total_tuples, sizeof(GOMO_RESULT_T), gomo_result_cmp); 
  //copy to ARRAY_T for processing by compute_qvalues
  ARRAY_T *score_array = allocate_array(bins->total_tuples);
  for (i = 0; i < bins->total_tuples; ++i) {
    set_array_item(i, results[i].empirical_pvalue, score_array);
  }
  // compute the qvalues
  compute_qvalues(
    false, // Don't stop with FDR.
    true,  // Do estimate pi-zero.
    NULL,  // Don't store pi-zero in a file.
    NUM_BOOTSTRAPS,
    NUM_BOOTSTRAP_SAMPLES,
    NUM_LAMBDA,
    MAX_LAMBDA,
    bins->total_tuples,
    score_array,
    NULL // No sampled p-values provided
  );
  //output motif's results to xml
  if (!text_only) {
    fprintf(xml_output,"\t<motif id=\"%s\" genecount=\"%d\"", xmlify(motif_id, buf, true), seq_scores->length);
    if (logo_path) fprintf(xml_output," logo=\"%s\"", logo_path);
    fprintf(xml_output,">\n");
  }
  // output the header line to the TSV output
  for (i = 0; i < bins->total_tuples; ++i) {
    SEQ_LIST_T *go_term_genes = go_map+(results[i].term_load_num);
    char *term_id = results[i].term_id;
    double score = results[i].gomo_score;
    double pvalue = results[i].empirical_pvalue;
    double qvalue = get_array_item(i, score_array);
    int annotated = go_term_genes->length;

    if (qvalue > q_threshold) break;

    fprintf(txt_output, "%s\t%s\t%1.3e\t%1.3e\t%1.3e\n",
	motif_id, term_id, score, pvalue, qvalue);
    if (! text_only) {
      fprintf(xml_output,"\t\t<goterm id=\"%s\" score=\"%1.3e\" "
          "pvalue=\"%1.3e\" qvalue=\"%1.3e\" annotated=\"%d\" group=\"\" "
          "nabove=\"\" nbelow=\"\" implied=\"u\"\n"
          "\t\t\t\tname=\"\"",
          xmlify(term_id, buf, true), score, pvalue, qvalue, annotated);
      if (best_genes_max > 0) {
        int genes_found = 0;
        int target = min(best_genes_max, annotated);
        fprintf(xml_output, ">\n");
        //sort the go map for this gene
        qsort(go_term_genes->seq_load_nums, go_term_genes->length, sizeof(int), intcmp);
        //scan down the list for the primary species and find the best genes until we reach max
        for (j = 0; j < seq_scores->length; ++j) {
          AMA_SCORE_T *score = seq_scores->list+j;
          if (binary_search(&(score->seq_load_num), go_term_genes->seq_load_nums, go_term_genes->length, sizeof(int), intcmp) >= 0) {
            fprintf(xml_output, "\t\t\t<gene id=\"%s\" rank=\"%d\" />\n",
                xmlify(score->seq_id, buf, true), j+1);
            if (++genes_found >= target) break;
          }
        }
        fprintf(xml_output, "\t\t</goterm>\n");
      } else {
        fprintf(xml_output, "/>\n");
      }
    }
  }
  if (! text_only) fprintf(xml_output,"\t</motif>\n");

  //clean up
  free_array(score_array);
  free(results);
  str_destroy(buf, false);
}

/*
 * Attempt to parse an integer from the string.
 * If the conversion fails then display the passed error message and exit.
 */
int parse_integer(char *str, char *no_conversion_msg, char *partial_conversion_msg, char *underflow_msg, char *overflow_msg) {
  char *end;
  long l;
  errno = 0;
  l = strtol(str, &end, 10);
  if (end == str) {
    if (no_conversion_msg) die(no_conversion_msg, str);
  } else if (*end != '\0') {
    if (partial_conversion_msg) die(partial_conversion_msg, end);
  }
  if (errno == ERANGE || l < INT_MIN || l > INT_MAX) {
    if (l < INT_MIN) { // underflow
      if (underflow_msg) die(underflow_msg, str);
    } else {
      if (overflow_msg) die(overflow_msg, str);
    }
  }
  return (int)l;
}

/*
 * Attempt to parse a double from the string.
 * If the conversion fails then display the passed error message and exit.
 */
double parse_double(char *str, char *no_conversion_msg, char *partial_conversion_msg, char *underflow_msg, char *overflow_msg) {
  char *end;
  double d;
  errno = 0;
  d = strtod(str, &end);
  if (end == str) { // no conversion
    if (no_conversion_msg) die(no_conversion_msg, str);
  } else if (*end != '\0') {
    if (partial_conversion_msg) die(partial_conversion_msg, end);
  }
  if (errno == ERANGE) {
    if (d == 0) { // underflow
      if (underflow_msg) die(underflow_msg, str);
    } else {
      if (overflow_msg) die(overflow_msg, str);
    }
  }
  return d;
}

/*****************************************************************************
 * process_cmd_line
 * Extract the required information out of the command line.
 *****************************************************************************/
void process_cmd_line(
    int argc, 
    char *argv[], 
    uint32_t *seed,
    char **out_dir, 
    bool *clobber, 
    bool *text_only,
    ARRAYLST_T **filter_motifs, 
    bool *use_E_values, 
    double *score_E_thresh,
    int *shuffle_scores,
    int *min_gene_count,
    double *q_threshold,
    char ** map_path,
    char ** dag_path,
    char ** motifs_path,
    ARRAYLST_T **ama_paths
  ) 
{
  // set defaults
  *out_dir = default_out_dir;
  *clobber = true;
  *text_only = false;
  *filter_motifs = NULL;
  *use_E_values = true;
  *score_E_thresh = 0.0;
  *shuffle_scores = DEFAULT_SHUFFLES;
  *min_gene_count = DEFAULT_MIN_GENE_COUNT;
  *q_threshold = DEFAULT_Q_VALUE_THRESHOLD;
  *dag_path = NULL;
  *motifs_path = NULL;

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  const int num_options = 15;
  cmdoption const motif_scan_options[] = {
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"dag", REQUIRED_VALUE},
    {"motifs", REQUIRED_VALUE},
    {"text", NO_VALUE},
    {"shuffle_scores", REQUIRED_VALUE},
    {"motif",REQUIRED_VALUE},
    {"score_E_thresh", REQUIRED_VALUE},
    {"t", REQUIRED_VALUE},
    {"min_gene_count", REQUIRED_VALUE},
    {"gs", NO_VALUE},
    {"seed", REQUIRED_VALUE},
    {"nostatus", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"version", NO_VALUE}
  };

  int option_index = 0;

  // Define the usage message.
  char usage[] =
    "Usage: gomo [options] <GO-term database> <scoring file>+\n"
    "\n"
    "   Options:\n"
    "     --o <output dir>      name of the directory for output. Will not replace\n" 
    "                             an existing directory; default: 'gomo_out'\n"
    "     --oc <output dir>     name of the directory for output. Will replace an\n"
    "                             existing directory; default: 'gomo_out'\n"
    "     --text                output tab separated values format to standard out\n"
    "                             and don't create output directory or files;\n"
    "                             default: create HTML and XML in <output dir>;\n"
    "     --dag <godag>         path to the optional Gene Ontology DAG file to be\n"
    "                             used for highlighting the specific terms in the\n"
    "                             GOMo HTML output; default: no highlighting\n"
    "     --motifs <motifs>     path to the optional motif file in MEME format\n"
    "                             used to generate (all of the) scoring file(s);\n"
    "                             used for adding sequence logos to HTML output;\n"
    "                             default: no logos in output;\n"
    "     --motif <id>          limit results to specified motif; option may \n"
    "                             be repeated; default: use all motifs\n"
    "     --shuffle_scores <n>  generate empirical null by shuffling the sequence\n"
    "                             to score assignments and computing statistics <n>\n"
    "                             times; default: <n>=1000\n"
    "     --t <q>               the q-value threshold considered significant;\n"
    "                             default: <q>=0.05, q >= 1.0 shows all results\n"
    "     --score_E_thresh <E>  the score E-value threshold above which the same\n"
    "                             rank is given to all sequences; \n"
    "                             default: no threshold\n"
    "     --min_gene_count <n>  only consider GO terms annotated with a at least\n"
    "                             <n> genes; default: <n> = 1).\n"
    "     --gs                  extract gene scores rather than p-values from the\n"
    "                             CisML input file(s); default: use CisML p-values\n"
    "     --seed <seed>         seed the random number generator; different\n"
    "                             values of <seed> will give slightly different\n"
    "                             outputs; default: <seed> is chosen randomly\n"
    "     --nostatus            don't print progress messages to stderr\n"
    "     --verbosity [1|2|3|4] control level of progress messages;\n"
    "                             1 = Quiet, 2 = Normal, 3 = High, 4 = Very High;\n"
    "                             default: 2\n"
    "     --version              print the version and exit\n"
    "\n"
    "   See doc/gomo.html for information on the formats of <GO-term database>\n"
    "   and <scoring file>.\n"
    "\n";

  // Parse the command line.
  if (simple_setopt(argc, argv, num_options, motif_scan_options) != NO_ERROR) {
    die("Error processing command line options: option name too long.\n");
  }

  while (true) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }
    //output directory if not already
    if (strcmp(option_name, "o") == 0) {
      *out_dir = option_value;
      *clobber = false;
    }
    //output to directory always
    else if (strcmp(option_name, "oc") == 0) {
      *out_dir = option_value;
      *clobber = true;
    }
    // dag
    else if (strcmp(option_name, "dag") == 0) {
      *dag_path = option_value;
    }
    // motif images
    else if (strcmp(option_name, "motifs") == 0) {
      *motifs_path = option_value;
    }
    else if (strcmp(option_name, "text") == 0) {
      *text_only = true;  
    }
    //only calculate for selected motifs
    else if (strcmp(option_name, "motif") == 0) {
      if (*filter_motifs == NULL) {
        *filter_motifs = arraylst_create();
      }
      arraylst_add(option_value, *filter_motifs);
    }
    //shuffle the ama sequence ids n times
    else if (strcmp(option_name, "shuffle_scores") == 0) {
      //test that the option is only digits
      int n = parse_integer(option_value, 
          "The option shuffle_scores was set to \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option shuffle_scores ended with \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option shuffle_scores was set to %s, it requires a "
              "number larger than zero.\n",
          "The option shuffle_scores was set to %s which "
              "could not be used because the value is too large to "
              "represent internally.\n");
      if (n < 1)
        die ("The option shuffle_scores was set to %d, it requires a "
            "number larger than zero.\n", n);
      *shuffle_scores = n;
    }
    //E-values (from ama) above the threshold are ranked the same
    else if (strcmp(option_name, "score_E_thresh") == 0) {
      double x = parse_double(option_value, 
          "The option score_E_thresh was set to \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option score_E_thresh ended with \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option score_E_thresh was set to %s, which "
              "could not be used because the value requires more precision "
              "then is possible to represent internally.\n",
          "The option score_E_thresh was set to %s which "
              "could not be used because the value is too large to "
              "represent internally.\n");
      if (x <= 0.0)
        die ("The option score-E-thresh was set to %f, it requires a "
            "number larger than zero.\n", x);
      *score_E_thresh = x;
    }
    //threshold of GOMo qvalues expected to report
    else if (strcmp(option_name, "t") == 0) {
      double x = parse_double(option_value, 
          "The option t was set to \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option t ended with \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option t was set to %s, which "
              "could not be used because the value requires more precision "
              "then is possible to represent internally.\n",
          "The option t was set to %s which "
              "could not be used because the value is too large to "
              "represent internally.\n");
      if (x < 0.0 || x > 1.0)
        die ("The option t was set to %f, it requires a "
            "number between 0.0 and 1.0.\n", x);
      *q_threshold = x;
    }
    //only calculate for go ids that are annotated with more genes than the minimum gene count 
    else if (strcmp(option_name, "min_gene_count") == 0) {
      int n = parse_integer(option_value, 
          "The option min_gene_count was set to \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option min_gene_count ended with \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option min_gene_count was set to %s, it requires a "
              "number larger than zero.\n",
          "The option min_gene_count was set to %s which "
              "could not be used because the value is too large to "
              "represent internally.\n");
      if (n < 1)
        die ("The option min_gene_count was set to %d, it requires a "
            "number larger than zero.\n", n);
      *min_gene_count = n;
    }
    //use gene scores instead of converting p-values into e-values
    else if (strcmp(option_name, "gs") == 0) {
      *use_E_values = false;
    }
    else if (strcmp(option_name, "seed") == 0) {
      if (sscanf(option_value, "%" SCNu32, seed) == 1) {
        mt_seed32new(*seed);
      }
    }
    //don't report any status information
    else if (strcmp(option_name, "nostatus") == 0) {
      status = false;
    }
    //the level of output
    else if (strcmp(option_name, "verbosity") == 0) {
      int n = parse_integer(option_value, 
          "The option verbosity was set to \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option verbosity ended with \"%s\" which could "
              "not be interpreted as a number.\n",
          "The option verbosity was set to %s, it requires a "
              "number between 1 and 4 inclusive.\n",
          "The option verbosity was set to %s it requires a "
              "number between 1 and 4 inclusive.\n");
      if (n < 1 || n > 5)
        die ("The option verbosity was set to %d, it requires a "
            "number between 1 and 4 inclusive.\n", n);
      verbosity = n; 
    }
    // print version and exit
    else if (strcmp(option_name, "version") == 0) {
      fprintf(stdout, VERSION "\n");
      exit(EXIT_SUCCESS);
    }
  }

  // if we had any motifs then sort them
  if (*filter_motifs != NULL) {
    arraylst_qsort(arraylst_compar_txt, *filter_motifs);
  }

  // Must have go-map directory and at least one scoring file
  if (argc < option_index + 2) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_FAILURE);
  }

  *map_path = argv[option_index++];

  if (!file_exists(*map_path)) {
    die("GO map file \"%s\" does not exist.\n", *map_path);
  }

  *ama_paths = arraylst_create();
  while (option_index < argc) { 
    char *score_file = argv[option_index++];
    if (!file_exists(score_file)) {
      die("Scoring file \"%s\" does not exist.\n", score_file);
    }
    arraylst_add(score_file,*ama_paths);
  }
}

/*
 * duplicate_file
 * copies the source file to the destination file
 */
bool duplicate_file(const char *source_file, const char *dest_file, bool warn) {
  BUF_T *buffer;
  FILE *fp_in, *fp_out;

  fp_in = fopen(source_file, "r");
  if (fp_in == NULL) {
    if (warn)
      fprintf(stderr, "Unable to open source file \"%s\" for duplicating to \"%s\", error given as %s\n", 
          source_file, dest_file, strerror(errno));
    return false;
  }

  fp_out = fopen(dest_file, "w");
  if (fp_out == NULL) {
    if (warn)
      fprintf(stderr, "Unable to open destination file \"%s\" for duplication of \"%s\", error given as %s\n", 
          dest_file, source_file, strerror(errno));
    return false;
  }

  buffer = buf_create(100);

  while (!feof(fp_in)) {
    buf_fread(buffer, fp_in);
    if (ferror(fp_in)) {
      if (warn)
        fprintf(stderr, "Unable to read data from source file \"%s\" for duplication to \"%s\", error given as %s\n", 
            source_file, dest_file, strerror(ferror(fp_in)));
      fclose(fp_out);
      fclose(fp_in);
      buf_destroy(buffer);
      return false;
    }
    buf_flip(buffer);
    buf_fwrite(buffer, fp_out);
    if (ferror(fp_out)) {
      if (warn)
        fprintf(stderr, "Unable to write data to destination file \"%s\" for duplication of \"%s\", error given as %s\n", 
            dest_file, source_file, strerror(ferror(fp_out)));
      fclose(fp_in);
      fclose(fp_out);
      buf_destroy(buffer);
      return false;
    }
    buf_compact(buffer);
  }
  if (fclose(fp_in)) {
    if (warn)
      fprintf(stderr, "Source file \"%s\" closed badly, error given as %s\n", source_file, strerror(errno));
    fclose(fp_out);
    buf_destroy(buffer);
    return false;
  }
  buf_flip(buffer);
  while (buf_remaining(buffer)) {
    buf_fwrite(buffer, fp_out);
    if (ferror(fp_out)) {
      if (warn)
        fprintf(stderr, "Unable to write data to destination file \"%s\" for duplication of \"%s\", error given as %s\n", 
            dest_file, source_file, strerror(ferror(fp_out)));
      fclose(fp_out);
      buf_destroy(buffer);
      return false;
    }
  }
  buf_destroy(buffer);
  if (fclose(fp_out)) {
    if (warn)
      fprintf(stderr, "Destination file \"%s\" closed badly, error given as %s\n", dest_file, strerror(errno));
    return false;
  }
  return true; 
}

/*************************************************************************
 * Entry point for gomo
 *************************************************************************/
int main(int argc, char *argv[]) 
{
  int i, j, shuffle_scores, min_gene_count, map_seq_count;
  double score_E_thresh, q_threshold;
  char *out_dir, *map_path, *txt_path, *xml_path, *res_dir, *dag_path, 
    *motifs_path, *gene_url, *motif_id;
  STR_T *logo_name, *logo_path;
  bool clobber, text_only, use_Evalues;
  RBTREE_T *go_ids, *seq_ids, *motif_scores;
  RBNODE_T *node;
  ARRAYLST_T *filter_motifs, *ama_paths;
  SEQ_LIST_T *map;
  SHUFFLER_T *shuffler;
  BINS_T *bins;
  AMA_GROUP_T *m_scores;
  FILE *txt_fp;
  FILE *xml_fp;
  clock_t c0, c1;
  uint32_t seed;
  bool output_err = false;
  char *commandline = NULL;

  //set seeds on pseudo random number generators
  seed = mt_seed();

  process_cmd_line(
      argc, 
      argv,
      &seed,
      &out_dir,
      &clobber,
      &text_only,
      &filter_motifs,
      &use_Evalues,
      &score_E_thresh,
      &shuffle_scores,
      &min_gene_count,
      &q_threshold,
      &map_path,
      &dag_path,
      &motifs_path,
      &ama_paths
    );

  DISPLAY_STATUS("Starting %s with random seed %" PRIu32 "\n", program_name, seed);

  // measure time
  c0 = clock();

  // load all ids for go terms and sequences
  map = load_map(map_path, min_gene_count, &gene_url, &go_ids, &seq_ids); 
  map_seq_count = rbtree_size(seq_ids);

  // load ama scores
  load_ama_scores(ama_paths, use_Evalues, score_E_thresh, filter_motifs, seq_ids, &motif_scores);

  // create shuffler
  shuffler = shuffler_create(map_seq_count, rbtree_size(seq_ids));
  
  if (!text_only) {
    // configure output
    if (create_output_directory(out_dir, clobber, (verbosity >= NORMAL_VERBOSE))) {
      // Failed to create output directory.
      exit(1);
    }
    //create the resource subdirectory
    res_dir = make_path_to_file(out_dir, RES_DIRNAME);
    if (create_output_directory(res_dir, clobber, (verbosity >= NORMAL_VERBOSE))) {
      // Failed to create resource directory
      exit(1);
    }
    // Create the name of the XML output file
    // "<dir>/XML_FILENAME" and open it for writing
    xml_path = make_path_to_file(out_dir, XML_FILENAME);
    xml_fp = fopen(xml_path, "w");
    // Create the name of the text output file
    // "<dir>/TSV_FILENAME" and open it for writing
    txt_path = make_path_to_file(out_dir, TSV_FILENAME);
    txt_fp = fopen(txt_path, "w");

    //start output
    commandline = start_xml_output(xml_fp, seed, argc, argv, map_path, ama_paths, 
        out_dir, clobber, text_only, use_Evalues, score_E_thresh, 
        min_gene_count, filter_motifs, shuffle_scores, q_threshold, 
        gene_url);
  } else {
    txt_fp = stdout;
    xml_path = NULL;
    xml_fp = NULL;
    res_dir = NULL;
  }

  // load motif PWMs
  if (motifs_path) {
    MREAD_T *mread;
    MOTIF_T *pwm;
    mread = mread_create(motifs_path, OPEN_MFILE, true);
    while (mread_has_motif(mread)) {
      pwm = mread_next_motif(mread);
      m_scores = (AMA_GROUP_T*)rbtree_get(motif_scores, get_motif_id(pwm));
      if (m_scores != NULL) {
        m_scores->motif = pwm;
      } else {
        destroy_motif(pwm);
      }
    }
    mread_destroy(mread);
  }

  // print the text header
  fprintf(txt_fp, "Motif_Identifier\tGO_Term_Identifier\tGOMo_Score\tp-value\tq-value\n");

  logo_name = str_create(10);
  logo_path = str_create(20);
  //Loop on motif_scores
  for (i = 1, node = rbtree_first(motif_scores); node != NULL; ++i, node = rbtree_next(node)) {
    motif_id = rbtree_key(node);
    m_scores = (AMA_GROUP_T*)rbtree_value(node);
    
    DISPLAY_STATUS("Processing motif %s # %d of %d \n", motif_id, i, rbtree_size(motif_scores));

    // create one all inclusive bin 
    bins = bins_create(0, NULL);
    
    shuffler_identity(shuffler);
    calculate_scores_and_bin(go_ids, map, m_scores, false, shuffler, bins);
    bins_sort(bins);

    DISPLAY_STATUS("Generating null model scores for motif %s \n", motif_id);
    
    for (j = 1; j <= shuffle_scores; ++j) {
      shuffler_shuffle(shuffler);
      calculate_scores_and_bin(go_ids, map, m_scores, true, shuffler, bins);
      SKIPABLE_STATUS("Generated %d%%\n", (int)(((double)j / (double)shuffle_scores) * 100.0));  
    }

    // create motif logo
    if (m_scores->motif) {
      DISPLAY_STATUS("Creating motif logo %s \n", motif_id);
      str_setf(logo_name, "logo%d", i);
      str_clear(logo_path);
      str_append_path(logo_path, 2, res_dir, str_internal(logo_name));
      CL_create1(m_scores->motif, false, false, "GOMo", str_internal(logo_path),
          false, true);
      str_append(logo_name, ".png", 4);
    }

    // calculate and output with qvalues
    calculate_qvalues_and_output(motif_id, str_internal(logo_name), bins, txt_fp, xml_fp,
        text_only, q_threshold, map, m_scores->species[0], 10);
    // clean up
    bins_destroy(bins);
  }

  // Finish the TSV output
  char *version_message = "# GOMo (Gene Ontology for Motifs): Version " VERSION " compiled on " __DATE__ " at " __TIME__ "\n";
  fprintf(txt_fp, "\n%s", version_message);
  fprintf(txt_fp, "# The format of this file is described at %s/%s.\n", SITE_URL, "doc/gomo-output-format.html");
  fprintf(txt_fp, "# %s\n", commandline);
  free(commandline);

  str_destroy(logo_name, false);
  str_destroy(logo_path, false);

  if (!text_only) {
    // finish output
    fprintf(xml_fp,"</gomo>\n");
    fclose(xml_fp);
    fclose(txt_fp);
  }

  DISPLAY_STATUS("Finished %s \n", program_name);

  /***************************************************************************
   * clean up
   **************************************************************************/
  // destroy go map
  for (i = 0; i < rbtree_size(go_ids); ++i) {
      free(map[i].seq_load_nums);
  }
  free(map);
  // destroy motif_scores
  rbtree_destroy(motif_scores);
  //destroy shuffler
  shuffler_destroy(shuffler);
  //destroy go terms and sequences
  rbtree_destroy(go_ids);
  rbtree_destroy(seq_ids);
  //destroy file lists
  arraylst_destroy(NULL, ama_paths);
  
  if (filter_motifs != NULL) arraylst_destroy(NULL, filter_motifs);

  /***************************************************************************
   *  post processing of xml
   **************************************************************************/
  if (!text_only) {
    if (dag_path) {
      GODAG_T *godag = load_go_dag(dag_path);
      rewrite_gomo_xml(out_dir, xml_path, godag); 
      destroy_go_dag(godag);
    }
    // html files become to big for non-constrainted q-values
    if (q_threshold <= 0.5) {
      char *stylesheet_path, *html_path;
      stylesheet_path = make_path_to_file(get_meme_data_dir(), HTML_STYLESHEET);
      html_path = make_path_to_file(out_dir, HTML_FILENAME);
      if (file_exists(stylesheet_path)) {
        if (! print_xml_filename_to_filename_using_stylesheet(
            xml_path,         /* path to XML input file IN */
            stylesheet_path,  /* path to MEME XSL stylesheet IN */
            html_path       /* path to HTML output file IN */
            ) 
          ) output_err = true;
      } else {
        output_err = true;
        if (verbosity >= NORMAL_VERBOSE) 
          fprintf(stderr, "Warning: could not find the stylesheet \"%s\" required for transformation of xml to html. Have you installed %s correctly?\n", 
            stylesheet_path, program_name);
      }
      free(stylesheet_path);
      free(html_path);
    } else {
      output_err = true;
      if (verbosity >= NORMAL_VERBOSE) 
        fprintf(stderr, "Warning: Reporting threshold for q-value (%g) is too permissive to allow HTML output due to size and memory constraints.\n",
            q_threshold);
    }
    free(xml_path);
    free(res_dir);
  }

  // measure time
  c1 = clock();
  if (verbosity >= NORMAL_VERBOSE) { // starting time
    fprintf(stderr,"cycles (CPU);            %ld cycles\n", (long) c1);
    fprintf(stderr,"elapsed CPU time:        %f seconds\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  }

  return(output_err ? 1 : 0);
}
