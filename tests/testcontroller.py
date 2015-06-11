import unittest


from gluon.globals import Request

execfile("applications/mbi/controllers/default.py", globals())

class TestDefaultController(unittest.TestCase):
    def setUp(self):
        self.request = Request(globals())


    def test01_dictHasForm(self):
        responseDict = algorithm()
        self.assertIn('form', responseDict)


    def test02_returnValLenIsZero(self):
        responseDict = algorithm()
        self.assertEqual(len(responseDict['d']),0)

    def test03_dictHasReturnVal(self):
        responseDict = algorithm()
        self.assertIn('d', responseDict)


suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(TestDefaultController))
unittest.TextTestRunner(verbosity=2).run(suite)