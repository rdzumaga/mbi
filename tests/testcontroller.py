import unittest


from gluon.globals import Request

execfile("applications/mbi/controllers/default.py", globals())

class TestDefaultController(unittest.TestCase):
    def setUp(self):
        self.request = Request(globals())


    def test01_nw_dictHasForm(self):
        responseDict = algorithm()
        self.assertIn('form', responseDict)


    def test02_nw_returnValLenIsZero(self):
        responseDict = algorithm()
        self.assertEqual(len(responseDict['d']),0)

    def test03_nw_dictHasReturnVal(self):
        responseDict = algorithm()
        self.assertIn('d', responseDict)

    def test04_sw_dictHasForm(self):
        responseDict = smithWaterman()
        self.assertIn('form', responseDict)


    def test05_sw_returnValLenIsZero(self):
        responseDict = smithWaterman()
        self.assertEqual(len(responseDict['d']),0)

    def test06_sw_dictHasReturnVal(self):
        responseDict = smithWaterman()
        self.assertIn('d', responseDict)


    def test07_fasta_dictHasForm(self):
        responseDict = fasta()
        self.assertIn('form', responseDict)


    def test08_fasta_returnValLenIsZero(self):
        responseDict = fasta()
        self.assertEqual(len(responseDict['d']),0)

    def test09_fasta_dictHasReturnVal(self):
        responseDict = fasta()
        self.assertIn('d', responseDict)

    def test10_fasta_dbProperlyLoaded(self):
        responseDict = fasta()
        self.assertIn('db', responseDict)

    def test11_fasta_dbHasLength(self):
        responseDict = fasta()
        self.assertNotEqual(len(responseDict['db']),0)



suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(TestDefaultController))
unittest.TextTestRunner(verbosity=2).run(suite)