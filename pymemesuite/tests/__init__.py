from . import (
    test_doctest,
    test_common,
    test_fimo,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_common))
    suite.addTests(loader.loadTestsFromModule(test_fimo))
    return suite
