from libmeme.meme cimport MODEL, DATASET, S_POINT, PRIORS

cdef extern from "em.h" nogil:
    cdef void em(MODEL* model, DATASET* dataset, DATASET* neg_dataset, S_POINT* s_point)
    cdef void m_step(MODEL* model, DATASET* dataset, PRIORS* priors, double wnsitesw)
    cdef void e_step(MODEL* model, DATASET* dataset)
    cdef double tcm_e_step(MODEL* model, DATASET* dataset)
    cdef double like_e_step(MODEL* model, DATASET* dataset)
    cdef void discretize(MODEL* model, DATASET* dataset, DATASET* neg_dataset, double (*E_STEP)(MODEL*, DATASET*))
    cdef void set_z(MODEL* model, DATASET* dataset)
    cdef void convert_theta_to_log(MODEL* model, DATASET* dataset)
    cdef int estep_maxima_nsites(
        MODEL* model,
        DATASET* primary,
        DATASET* constol,
        double (*E_STEP)(MODEL*, DATASET*),
        bint primary_groups[3],
        bint control_groups[3],
        int min_nsites,
        int max_nsites,
        double thresh_in,
        double* log_pv,
        double* llr,
        double* thresh_out
    )
