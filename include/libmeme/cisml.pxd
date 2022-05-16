from libcpp cimport bool

from libc.stdio cimport FILE

from libmeme.array cimport ARRAY_T


cdef extern from "cisml.h" nogil:

    cdef struct cisml:
        pass
    ctypedef cisml CISML_T

    cdef struct cisml_match_it:
        pass
    ctypedef cisml_match_it CISML_MATCH_IT_T

    cdef struct multi_pattern:
        pass
    ctypedef multi_pattern MULTI_PATTERN_T

    cdef struct multi_pattern_match:
        pass
    ctypedef multi_pattern_match MULTI_PATTERN_MATCH_T

    cdef struct pattern:
        pass
    ctypedef pattern PATTERN_T

    cdef struct scanned_sequence:
        pass
    ctypedef scanned_sequence SCANNED_SEQUENCE_T

    cdef struct matched_element:
        pass
    ctypedef matched_element MATCHED_ELEMENT_T

    cdef CISML_T* allocate_cisml(
        const char* program_name,
        const char* command_line,
        const char* pattern_file,
        const char* sequence_file
    )
    cdef void free_cisml(CISML_T* cisml)

    cdef char *get_cisml_program_name(CISML_T *cisml)
    cdef char *get_cisml_command_line(CISML_T *cisml)
    cdef char *get_cisml_pattern_file(CISML_T *cisml)
    cdef char *get_cisml_sequence_file(CISML_T *cisml)

    cdef void set_cisml_background_file(CISML_T *cisml, char *background_file)
    cdef char *get_cisml_background_file(CISML_T *cisml)

    cdef void set_cisml_pattern_pvalue_cutoff(CISML_T *cisml, double pattern_pvalue_cutoff)
    cdef void clear_cisml_pattern_pvalue_cutoff(CISML_T *cisml)
    cdef bool has_cisml_pattern_pvalue_cutoff(CISML_T *cisml)
    cdef double get_cisml_pattern_pvalue_cutoff(CISML_T *cisml)

    cdef void set_cisml_sequence_pvalue_cutoff(CISML_T *cisml, double sequence_pvalue_cutoff)
    cdef void clear_cisml_sequence_pvalue_cutoff(CISML_T *cisml)
    cdef bool has_cisml_sequence_pvalue_cutoff(CISML_T *cisml)
    cdef double get_cisml_sequence_pvalue_cutoff(CISML_T *cisml)

    cdef int get_cisml_num_passing_cutoff(CISML_T *cisml)
    cdef void set_cisml_site_pvalue_cutoff(CISML_T *cisml, double site_pvalue_cutoff)
    cdef void clear_cisml_site_pvalue_cutoff(CISML_T *cisml)
    cdef bool has_cisml_site_pvalue_cutoff(CISML_T *cisml)
    cdef double get_cisml_site_pvalue_cutoff(CISML_T *cisml)

    cdef void set_cisml_site_qvalue_cutoff(CISML_T *cisml, double site_qvalue_cutoff)
    cdef void clear_cisml_site_qvalue_cutoff(CISML_T *cisml)
    cdef bool has_cisml_site_qvalue_cutoff(CISML_T *cisml)
    cdef double get_cisml_site_qvalue_cutoff(CISML_T *cisml)

    cdef void set_cisml_sequence_filter(CISML_T *cisml, char *sequence_filter)
    cdef char *get_cisml_sequence_filter(CISML_T *cisml)

    cdef int get_cisml_num_multi_patterns(CISML_T *cisml)
    cdef MULTI_PATTERN_T **get_cisml_multi_patterns(CISML_T *cisml)

    cdef int get_cisml_num_patterns(CISML_T *cisml)
    cdef PATTERN_T **get_cisml_patterns(CISML_T *cisml)

    cdef int get_cisml_num_stored_matches(CISML_T *cisml)

    cdef void add_cisml_multi_pattern(CISML_T *cisml, MULTI_PATTERN_T* multi_pattern)
    cdef void add_cisml_pattern(CISML_T *cisml, PATTERN_T* pattern)

    cdef CISML_MATCH_IT_T *allocate_cisml_match_iterator(CISML_T *cisml)
    cdef void free_cisml_match_iterator(CISML_MATCH_IT_T *it)
    cdef MATCHED_ELEMENT_T *cisml_match_iterator_next(CISML_MATCH_IT_T *it)

    cdef MULTI_PATTERN_T *allocate_multi_pattern()
    cdef void free_multi_pattern(MULTI_PATTERN_T *multi_pattern)

    cdef void set_multi_pattern_score(MULTI_PATTERN_T *multi_pattern, double score)
    cdef bool has_multi_pattern_score(MULTI_PATTERN_T *multi_pattern)
    cdef void clear_multi_pattern_score(MULTI_PATTERN_T *multi_pattern)
    cdef double get_multi_pattern_score(MULTI_PATTERN_T *multi_pattern)

    cdef void set_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern, double pvalue)
    cdef bool has_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern)
    cdef void clear_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern)
    cdef double get_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern)

    cdef int get_multi_pattern_num_patterns(MULTI_PATTERN_T *multi_pattern)
    cdef MULTI_PATTERN_MATCH_T *get_multi_pattern_match(MULTI_PATTERN_T *multi_pattern)
    cdef void set_multi_pattern_match(MULTI_PATTERN_T *multi_pattern,  MULTI_PATTERN_MATCH_T *match)

    cdef PATTERN_T **get_multi_pattern_patterns(MULTI_PATTERN_T *multi_pattern)
    cdef void add_multi_pattern_pattern(MULTI_PATTERN_T *multi_pattern, PATTERN_T* pattern)
    cdef void multi_pattern_calculate_qvalues(
      int num_multi_pattern,
      MULTI_PATTERN_T **multi_patterns,
      ARRAY_T *sampled_pvalues
    )

    cdef MULTI_PATTERN_MATCH_T *allocate_multi_pattern_match(
        char *seq_name,
        char *seq,
        int start,
        int stop,
        double pvalue
    )
    cdef void free_multi_pattern_match(MULTI_PATTERN_MATCH_T *match)

    cdef void set_multi_pattern_match_evalue(MULTI_PATTERN_MATCH_T *match, double evalue)
    cdef double get_multi_pattern_match_pvalue(MULTI_PATTERN_MATCH_T* match)

    cdef void set_multi_pattern_match_qvalue(MULTI_PATTERN_MATCH_T *match, double qvalue)
    cdef double get_multi_pattern_match_qvalue(MULTI_PATTERN_MATCH_T* match)

    cdef PATTERN_T *allocate_pattern(char *accession, char *name)
    cdef void free_pattern(PATTERN_T *pattern)

    cdef void add_pattern_elements_to_scanned_seq(PATTERN_T *pattern)
    cdef bool add_pattern_matched_element(PATTERN_T *pattern, MATCHED_ELEMENT_T *element)
    cdef bool set_pattern_max_stored_matches(PATTERN_T *pattern, int max)
    cdef void set_pattern_max_pvalue_retained(PATTERN_T *pattern, double max_pvalue)
    cdef int get_pattern_max_stored_matches(PATTERN_T *pattern)
    cdef int get_pattern_num_stored_matches(PATTERN_T *pattern)
    cdef bool get_pattern_is_complete(PATTERN_T *pattern)
    cdef void set_pattern_is_complete(PATTERN_T *pattern)
    cdef double get_pattern_max_pvalue_retained(PATTERN_T *pattern)
    cdef char *get_pattern_name(PATTERN_T* pattern)
    cdef char *get_pattern_accession(PATTERN_T* pattern)

    cdef void set_pattern_pvalue(PATTERN_T *pattern, double pvalue)
    cdef void clear_pattern_pvalue(PATTERN_T *pattern)
    cdef bool has_pattern_pvalue(PATTERN_T *pattern)
    cdef double get_pattern_pvalue(PATTERN_T* pattern)

    cdef bool has_pattern_qvalues(PATTERN_T *pattern)
    cdef void set_pattern_score(PATTERN_T *pattern, double score)
    cdef void clear_pattern_score(PATTERN_T *pattern)
    cdef bool has_pattern_score(PATTERN_T *pattern)
    cdef double get_pattern_score(PATTERN_T* pattern)
    cdef void set_pattern_db(PATTERN_T *pattern, char *db)
    cdef char *get_pattern_db(PATTERN_T* pattern)
    cdef void set_pattern_lsid(PATTERN_T *pattern, char *lsid)
    char *get_pattern_lsid(PATTERN_T* pattern)

    cdef int get_pattern_num_scanned_sequences(PATTERN_T *pattern)
    cdef long get_pattern_num_scanned_positions(PATTERN_T *pattern)
    cdef int get_pattern_num_stored_matches(PATTERN_T *pattern)

    cdef MATCHED_ELEMENT_T **get_pattern_stored_matches(PATTERN_T *pattern)
    cdef MATCHED_ELEMENT_T **get_pattern_matched_elements(PATTERN_T *pattern)

    cdef void pattern_calculate_qvalues(PATTERN_T *pattern, ARRAY_T *sampled_pvalues)
    cdef SCANNED_SEQUENCE_T **get_pattern_scanned_sequences(PATTERN_T *pattern)
    cdef SCANNED_SEQUENCE_T *allocate_scanned_sequence(char *accession, char *name, PATTERN_T *parent)
    cdef void free_scanned_sequence(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef void free_scanned_sequence_from_matched_elements(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef PATTERN_T *get_scanned_sequence_parent(SCANNED_SEQUENCE_T *scanned_sequence)

    cdef char *get_scanned_sequence_name(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef char *get_scanned_sequence_accession(SCANNED_SEQUENCE_T* scanned_sequence)

    cdef void set_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence, double pvalue)
    cdef void clear_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef bool has_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef double get_scanned_sequence_pvalue(SCANNED_SEQUENCE_T* scanned_sequence)

    cdef void set_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence, double score)
    cdef void clear_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef bool has_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef double get_scanned_sequence_score(SCANNED_SEQUENCE_T* scanned_sequence)

    cdef void set_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence, int length)
    cdef void clear_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef bool has_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence)
    cdef int get_scanned_sequence_length(SCANNED_SEQUENCE_T* scanned_sequence)

    cdef void set_scanned_sequence_db(SCANNED_SEQUENCE_T *scanned_sequence, char *db)
    cdef char *get_scanned_sequence_db(SCANNED_SEQUENCE_T* scanned_sequence)

    cdef void set_scanned_sequence_lsid(SCANNED_SEQUENCE_T *scanned_sequence, char *lsid)
    cdef char *get_scanned_sequence_lsid(SCANNED_SEQUENCE_T* scanned_sequence)

    cdef void add_scanned_sequence_matched_element(SCANNED_SEQUENCE_T *sequence, MATCHED_ELEMENT_T *element)
    cdef void add_scanned_sequence_scanned_position(SCANNED_SEQUENCE_T *sequence)
    cdef int get_scanned_sequence_num_matched_elements(SCANNED_SEQUENCE_T *sequence)
    cdef long get_scanned_sequence_num_scanned_positions(SCANNED_SEQUENCE_T *sequence)

    cdef int get_scanned_sequence_num_matched_elements(SCANNED_SEQUENCE_T *sequence)
    cdef long get_scanned_sequence_num_scanned_positions(SCANNED_SEQUENCE_T *sequence)

    cdef MATCHED_ELEMENT_T **get_scanned_sequence_matched_elements(SCANNED_SEQUENCE_T *sequence)
    cdef MATCHED_ELEMENT_T **get_scanned_sequence_matched_elements(SCANNED_SEQUENCE_T *sequence)
    cdef MATCHED_ELEMENT_T *allocate_matched_element(
        int start,
        int stop,
        SCANNED_SEQUENCE_T *parent
    )
    cdef MATCHED_ELEMENT_T *allocate_matched_element_without_inversion(
        int start,
        int stop,
        const char *seq,
        SCANNED_SEQUENCE_T *parent
    )
    cdef MATCHED_ELEMENT_T *allocate_matched_element_with_score(
        int start,
        int stop,
        double score,
        double pvalue,
        SCANNED_SEQUENCE_T *parent
    )
    cdef void free_matched_element(MATCHED_ELEMENT_T *element)


    cdef int get_matched_element_start(MATCHED_ELEMENT_T* matched_element)
    cdef void set_matched_element_start(MATCHED_ELEMENT_T* matched_element, int newstart)

    cdef int get_matched_element_stop(MATCHED_ELEMENT_T* matched_element)
    cdef void set_matched_element_stop(MATCHED_ELEMENT_T* element, int newstop)

    cdef void set_matched_element_score(MATCHED_ELEMENT_T *element, double score)
    cdef bool has_matched_element_score(MATCHED_ELEMENT_T *element)
    cdef double get_matched_element_score(MATCHED_ELEMENT_T* element)

    cdef void set_matched_element_pvalue(MATCHED_ELEMENT_T *element, double pvalue)
    cdef bool has_matched_element_pvalue(MATCHED_ELEMENT_T *element)
    cdef double get_matched_element_pvalue(MATCHED_ELEMENT_T* element)

    cdef bool has_matched_element_qvalue(MATCHED_ELEMENT_T *element)
    cdef double get_matched_element_qvalue(MATCHED_ELEMENT_T* element)

    cdef void set_matched_element_clusterid(MATCHED_ELEMENT_T *element, char *clusterid)
    cdef char *get_matched_element_clusterid(MATCHED_ELEMENT_T* element)

    cdef const char *get_matched_element_sequence(MATCHED_ELEMENT_T* element)

    cdef void set_matched_element_strand(MATCHED_ELEMENT_T* element, char strand)
    cdef char get_matched_element_strand(MATCHED_ELEMENT_T* element)

    cdef SCANNED_SEQUENCE_T *get_matched_element_scanned_seq(MATCHED_ELEMENT_T* element)
    cdef void set_matched_element_sequence(MATCHED_ELEMENT_T* element, char *seq)

    cdef void sort_matched_elements(
        bool sort_by_pvalue,
        int num_elements,
        MATCHED_ELEMENT_T **elements
    )

    cdef void print_cisml_start(
        FILE* out,
        char *program_name,
        bool print_header,
        const char *stylesheet,
        bool print_namespace
    )
    cdef void print_cisml_end(FILE* out)
    cdef void print_cisml_parameters(FILE *out, CISML_T *cisml)
    cdef void print_cisml(
    		FILE* out,
    		CISML_T *cisml,
    		bool print_header,
    		const char *stylesheet,
    		bool print_namespace
    )
    cdef void print_cisml_start_pattern(
        CISML_T *cisml,
        FILE *out,
        PATTERN_T *pattern
    )
    cdef void print_cisml_end_pattern(FILE *out)
    cdef void print_cisml_scanned_sequences(
        CISML_T *cisml,
        FILE *out,
        int num_seqs,
        SCANNED_SEQUENCE_T **sequences
    )
    cdef void print_cisml_scanned_sequence_start(
        CISML_T *cisml,
        FILE *out,
        SCANNED_SEQUENCE_T *seq
    )
    cdef void print_cisml_scanned_sequence_end(FILE *out)

    cdef bool print_full_results(
        CISML_T *cisml,
        char *output_dirname,
        char *xml_filename,
        char *html_filename,
        char *text_filename,
        char *gff_filename,
        bool allow_clobber,
        bool print_namespace
    )
    cdef bool print_cisml_as_text(CISML_T *cisml)
    cdef CISML_T* read_cisml(const char* cisml_filename)
