from libc.stdint cimport uint8_t

from libmeme.alphabet cimport ALPH_T, ALPHABET_T
from libmeme.array cimport ARRAY_T
from libmeme.data_types cimport Z_T, LCB_T, WEIGHTS_T
from libmeme.hash_table cimport HASH_TABLE
from libmeme.heap cimport HEAP
from libmeme.mtype cimport MOTYPE
from libmeme.prior cimport PriorLib

cdef extern from "user.h" nogil:

    const int MAXALPH
    const int MAXSITE


cdef extern from "meme.h" nogil:

    ctypedef enum MASK_TYPE:
        Sym
        Tri

    ctypedef enum MASK_COMB_TYPE:
        Minp
        Minpop

    ctypedef enum NEGTYPE:
        Pair
        Blend

    ctypedef enum MAP_TYPE:
        Uni
        Pam

    ctypedef enum PTYPE:
        Mega
        MegaP
        Dmix
        Dirichlet
        Addone

    ctypedef enum STYPE:
        Combine
        Separate
        Norc
        Unstranded

    ctypedef enum OBJTYPE:
        Classic
        NC
        SE
        DE
        NZ
        CE
        CD

    ctypedef enum TESTTYPE:
        mHG
        mRS
        mBN
        No_testtype

    ctypedef enum GLOB_OBJ_FN:
        LLR_POP
        LLR
        Likelihood
        Class
        GLAM

    ctypedef enum POINT_BRANCHES:
        NORMAL
        X_ONLY
        ALL
        NO_POINT_B

    const size_t MAX_LOCAL_MAXIMA

    cdef struct sample:
        char*      sample_name
        char*      descript
        long       length
        char*      seq
        uint8_t*   res
        uint8_t*   resic
        double     sw
        WEIGHTS_T* weights
        double*    not_o
        int*       log_not_o
        int**      pY
        char*      pYic
        Z_T*       z
        Z_T*       z_buf
        double*    psp_original
        double*    psp_original_buf
        double*    log_psp
        double     max_log_psp
        double*    log_psp_buf
        double*    counts
        LCB_T*     logcumback
        int        nsites
        int*       sites
        double     minpv
        int        group
        int        orig_index
    ctypedef sample SAMPLE

    ctypedef double** THETA

    ctypedef struct p_prob:
        int    x
        int    y
        bint   ic
        bint   negative
        int    rank
        double ranksum
        double dist
        double in_region
        double wN
        double prob
    ctypedef p_prob* P_PROB

    cdef struct Model:
        MOTYPE  mtype
        int     min_w
        int     max_w
        bint    all_widths
        double  pw
        double  psites
        P_PROB  maxima
        bint    pal
        bint    invcomp
        int     imotif
        int     w
        THETA   theta
        THETA   logtheta
        THETA   logtheta_rc
        THETA   obs
        double  lambda_      "lambda"
        double  lambda_obs
        double  nsites
        double  nsites_obs
        int     nsites_dis
        char    cons[MAXSITE+1]
        char    cons0[MAXSITE+1]
        double  rentropy[MAXSITE+1]
        double  rel
        double  ic
        double  ll
        double  mll_0
        double  mll_1
        double  logpv
        double  logev
        double  llr
        double  site_threshold
        int     iter
        int     ID
        int     iseq
        int     ioff
        double  pc
        OBJTYPE objfun
        int     alength
    ctypedef Model MODEL

    cdef struct p_point:
        int       c
        int*      w
        double*   nsites
        uint8_t** e_cons0
    ctypedef p_point P_POINT

    cdef struct s_point:
        double   score
        int      iseq
        int      ioff
        int      w0
        double   nsites0
        double   wgt_nsites
        uint8_t* e_cons0
        char*    cons0
        HEAP*    seed_heap
        bint     evaluate
        double   sig
    ctypedef s_point S_POINT

    cdef struct candidate:
        S_POINT* s_point
        int      w
        bint     pal
        bint     invcomp
        double   lambda_ "lambda"
        char     cons[MAXSITE+1]
        double   ic
        double   rel
        double   ll
        double   sig
        double   llr
    ctypedef candidate CANDIDATE

    cdef struct motif_summary:
        int      width
        int      num_sites
        int      num_negative_sites
        double   ic
        double   re
        double   llr
        double   p_value_exp
        double   p_value_mant
        double   e_value_exp
        double   e_value_mant
        double   bayes
        double   elapsed_time
        double** pssm
        THETA    psfm
        char*    regexp
        char*    consensus
        P_PROB   sites
    ctypedef motif_summary MOTIF_SUMMARY

    cdef struct Priors:
        PTYPE     ptype
        double    prior_count[MAXALPH]
        PriorLib* plib
        PriorLib* plib0
    ctypedef Priors PRIORS

    cdef struct motif:
        char*      name
        int        width
        int        pos
        double     roc
        int        shift
        int        pass_ "pass"
        double     recall
        double     precision
        double     min_thresh
        double     max_thresh
        double     pal
        double     like
        double     sig
        double     ic
        double     sites
        int        w
        double     thresh
        HASH_TABLE ht
    ctypedef motif MOTIF

    cdef struct branch_params:
        int            bfactor
        POINT_BRANCHES point_branch
        bint           w_branch
    ctypedef branch_params BRANCH_PARAMS

    cdef struct hist_itm:
        int x
        int count
    ctypedef hist_itm HIST_ITM

    cdef struct hist:
        int       n
        HIST_ITM* entries
    ctypedef hist HIST

    cdef struct Group_t:
        bint em[3]
        bint trim[3]
        bint pvalue[3]
        bint nsites[3]
    ctypedef Group_t GROUP_T

    cdef struct Dataset:
        ALPH_T*  alph
        int      total_res
        double   wgt_total_res
        int      n_samples
        int      n_group[3]
        int      group_last_idx[3]
        SAMPLE** samples
        SAMPLE** input_order
        double*  seq_weights
        int      n_wgts
        long     max_slength
        long     min_slength
        double   ce_frac
        int      n_region[3]
        int      region_last_pos[3]
        int      ce_max_dist
        int      psp_w
        int      log_psp_w
        double*  res_freq
        bint     mpi
        bint     invcomp
        bint     pal
        THETA    map
        THETA    lomap
        MOTIF*   motifs
        NEGTYPE  negtype
        int      back_order
        ARRAY_T* back
        double   log_total_prob
        PRIORS*  priors
        P_POINT* p_point
        double   min_nsites
        double   max_nsites
        double   wnsites
        bint     ma_adj
        double   wg
        double   ws
        bint     endgaps
        double   distance
        int      nmotifs
        int      maxiter
        double   evt
        char*    mod
        char*    mapname
        double   map_scale
        char*    priorname
        double   beta
        int      seed
        double   seqfrac
        char*    plib_name
        char*    bfile
        char*    datafile
        char*    negfile
        int      neg_n_samples
        char*    pspfile
        char*    command
        OBJTYPE  objfun
        THETA    pairwise
        char*    output_directory
        double   max_time
        int      main_hs
        TESTTYPE test
        double   hsfrac
        int      shuffle
        double   hs_decrease
        BRANCH_PARAMS* branch_params
        bint     use_llr
        int      brief
        bint     print_heaps
        bint     print_pred
        bint     include_group[3]
        bint     save_include_group[3]
        bint     include_region[3]
        bint     save_include_region[3]
        int      n_included
        GROUP_T  primary_groups
        GROUP_T  control_groups
        int      last_seed_seq
        int      last_seed_pos
        int      search_size
        int      max_words
        int      max_size
        int      no_rand
        int      classic_max_nsites
        int      imotif
    ctypedef Dataset DATASET

    ctypedef struct SITE:
        int    seqno
        int    pos
        double zij
        int    nvcomp

    ctypedef struct TILING:
        int*    hits
        double* pvalues
        int*    svalues
        double  pv
        char*   diagram

    cdef void set_seq_groups_to_include(DATASET* dataset, bint groups[3])
    cdef void set_seq_regions_to_include(DATASET* dataset, bint incl0, bint incl1, bint incl2)

    cdef void copy_model(MODEL *m1, MODEL *m2, ALPH_T *alph)
    cdef void free_model(MODEL *model)
    cdef MODEL* create_model(MOTYPE mtype, bint invcomp, int max_w, ALPH_T* alph, OBJTYPE objfun)

    cdef S_POINT *get_starts(DATASET *primary, DATASET* control, MODEL* model, uint8_t* e_cons, int* n_starts)
    cdef THETA init_map(MAP_TYPE type, double scale, ALPH_T* alph, ARRAY_T* back, bint lo)
    cdef void convert_to_lmap(THETA map, int lmap[MAXALPH][MAXALPH], ALPH_T* alph)

    cdef bint init_model(S_POINT* s_point, MODEL* model, DATASET* dataset, int imotif)
    cdef void erase(DATASET* dataset, MODEL* model)
    cdef PRIORS *create_priors(PTYPE ptype, double beta, DATASET* dataset, char* plib_name)
    cdef void init_meme_background(char* bfile, bint rc, DATASET* dataset, char* alph_file, ALPHABET_T alphabet_type, int order, char* seqfile, bint status)
