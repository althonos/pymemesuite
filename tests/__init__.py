from . import (
    test_doctest,
    test_common,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_common))
    return suite
