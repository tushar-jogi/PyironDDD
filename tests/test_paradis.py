from PyironDDD.paradisjob import ParaDis
import unittest


class TestParadis(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))


if __name__ == '__main__':
    unittest.main()
