import unittest
import parse_data
class TestModelClass(unittest.TestCase):

    def test0_read_msa_and_encoding(self):
        encoding_seq=parse_data.read_msa_and_encoding(file_name="test")
        self.assertEqual(len(encoding_seq.columns),6532)