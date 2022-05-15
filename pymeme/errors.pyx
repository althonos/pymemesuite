# coding: utf-8
# cython: language_level=3, linetrace=True
"""Common errors for the MEME modules.
"""

class AllocationError(MemoryError):
    """A memory error that is caused by an unsuccessful allocation.
    """

    def __init__(self, str ctype):
        self.ctype = ctype

    def __repr__(self):
        return "{}({!r})".format(type(self).__name__, self.ctype)

    def __str__(self):
        return "Could not allocate {!r}".format(self.ctype)
