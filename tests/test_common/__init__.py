from . import (
    test_alphabet,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_alphabet))
    return suite
