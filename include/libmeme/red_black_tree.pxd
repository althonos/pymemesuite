cdef extern from "red-black-tree.h" nogil:

    cdef struct rbtree_t:
        pass
    ctypedef rbtree_t RBTREE_T

    cdef struct rbnode_t:
        pass
    ctypedef rbnode_t RBNODE_T

    cdef void rbtree_check(RBTREE_T* tree)
