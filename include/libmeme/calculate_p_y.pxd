from libmeme.alphabet cimport ALPH_T
from libmeme.meme cimport DATASET, MAXALPH
from libmeme.sp_matrix cimport SP_MATRIX
from libmeme.user cimport MAXSITE, DATASET, MODEL


cdef extern from "calculate_p_y.h" nogil:
  cdef void get_pY(DATASET* dataset, int* theta_1[MAXSITE], int w, int pYindex)
  cdef void combine_strands(DATASET* dataset, int w)

  cdef void create_lmotif(ALPH_T* alph, char* seed_str, int lmap[MAXALPH][MAXALPH], int* lmotif, int* moti_width)
  cdef void init_lmotif(int w, uint8_t res, int* theta_1[MAXSITE], int lmap[MAXALPH][MAXALPH])

  cdef void evaluate_seed_DP(char* new_seed, SEED_DIFFS* s_diffs, int lmap[MAXALPH][MAXALPH], MOTYPE mtype, bint ic, DATASET* dataset, SP_MATRIX* sp_matrix)
  cdef void next_pY_branching(int* lmotif_new[MAXSITE], SEED_DIFFS* s_diffs, DATASET* dataset, int pYindex)
