import unittest
from sequence import Sequences


class TestModelClass(unittest.TestCase):

    seq = Sequences()
    def test0_get_AA_within_distance_from_structure_file(self):
        """
        This is the test function for get_AA_within_distance_from_structure_file() in Sequences class

        """
        acs_dictionary, seq_length = self.seq.get_AA_within_distance_from_structure_file(
            residue_dictionary={11:"Y"})
        print(acs_dictionary)
        list1 = acs_dictionary[11]
        self.assertEqual(9, list1[0][0])
        self.assertEqual('ASP', list1[0][1])
        self.assertEqual(int(4.878653-list1[0][2]), 0)
        self.assertEqual(seq_length,261)

    def test1_amino_acid_properties(self):
        """
        This is the function for test amino_acid_properties() in Sequences class
        :return:
        """
        amino_acid = "Y"
        properties = self.seq.amino_acid_properties(amino_acid, "Y")
        self.assertEqual(properties["charge"],0)
        self.assertEqual(properties["volume"]-6.47, 0)
        self.assertEqual(properties["similarity_score"] - 7.0, 0)

    def test2_remove_sequences_from_result_of_mmseqs(self):
        self.seq.remove_sequences_from_result_of_mmseqs(
            "../autodata/sequences/alnRes.tab",
            seq_file="../autodata/sequences/PF08241.15PF03602.18.fasta")

def main():
    unittest.main()
    print(0)


if __name__ == "__main__":
    main()