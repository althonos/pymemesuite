from libmeme.array cimport ARRAY_T

cdef extern from "prior-dist.h" nogil:

    cdef struct prior_dist:
        pass
    ctypedef prior_dist PRIOR_DIST_T

    PRIOR_DIST_T* new_prior_dist(const char* filename)
    void free_prior_dist(PRIOR_DIST_T* prior_dist)
    double get_prior_dist_minimum(PRIOR_DIST_T* prior_dist)
    double get_prior_dist_maximum(PRIOR_DIST_T* prior_dist)
    double get_prior_dist_median(PRIOR_DIST_T* prior_dist)
    ARRAY_T* get_prior_dist_array(PRIOR_DIST_T* prior_dist)
    int get_prior_dist_length(PRIOR_DIST_T* prior_dist)
    double get_prior_dist_offset(PRIOR_DIST_T* prior_dist)
    double get_prior_dist_scale(PRIOR_DIST_T* prior_dist)
    double get_min_lo_prior(PRIOR_DIST_T* prior_dist, double alpha)
    double get_max_lo_prior(PRIOR_DIST_T* prior_dist, double alpha)
