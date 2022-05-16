

cdef extern from "reservoir.h" nogil:

    cdef struct reservoir_sampler:
       pass
    ctypedef reservoir_sampler RESERVOIR_SAMPLER_T

    cdef RESERVOIR_SAMPLER_T* new_reservoir_sampler(size_t size, void(*free_sample)(void*))

    cdef void free_reservoir(RESERVOIR_SAMPLER_T* reservoir)
    cdef void clear_reservoir(RESERVOIR_SAMPLER_T* reservoir)

    cdef size_t get_reservoir_size(RESERVOIR_SAMPLER_T* reservoir)
    cdef size_t get_reservoir_num_samples_seen(RESERVOIR_SAMPLER_T* reservoir)
    cdef size_t get_reservoir_num_samples_retained(RESERVOIR_SAMPLER_T* reservoir)
    cdef double* get_reservoir_samples(RESERVOIR_SAMPLER_T* reservoir)
    cdef void reservoir_sample(RESERVOIR_SAMPLER_T* reservoir, double value)
    cdef void reservoir_sample_pointer(RESERVOIR_SAMPLER_T* reservoir, void* sample)
