from . import (
    test_alphabet,
    test_background,
    test_motif_file,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_alphabet))
    suite.addTests(loader.loadTestsFromModule(test_background))
    suite.addTests(loader.loadTestsFromModule(test_motif_file))
    return suite
