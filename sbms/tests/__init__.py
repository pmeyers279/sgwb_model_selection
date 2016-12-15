import unittest
def get_suite():
    import sbms.tests
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(sbms.tests)
    return suite
