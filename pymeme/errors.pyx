# coding: utf-8
# cython: language_level=3, linetrace=True
"""Common errors for the MEME modules.
"""

class AllocationError(MemoryError):
    """A memory error that is caused by an unsuccessful allocation.
    """

    def __init__(self, str ctype, int itemsize, int count=1):
        self.ctype = ctype
        self.itemsize = itemsize
        self.count = count

    def __repr__(self):
        cdef str typename = type(self).__name__
        if self.count == 1:
            return f"{typename}({self.ctype!r}, {self.itemsize})"
        return f"{typename}({self.ctype!r}, {self.itemsize}, {self.count})"

    def __str__(self):
        if self.count == 1:
            return f"Could not allocate {self.itemsize} bytes for type {self.ctype}"
        return f"Could not allocate {self.itemsize*self.count} bytes for an array of {self.count} {self.ctype}"
