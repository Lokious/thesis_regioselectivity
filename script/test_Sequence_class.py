import unittest
from sequence import Sequences
class TestModel_class(unittest.TestCase):

    seq = Sequences()
    def test0_get_AA_within_distance_from_structure_file(self):
        """
        This is the test function for get_AA_within_distance_from_structure_file in Sequences class
    
        """
        acs_dictionary=self.seq.get_AA_within_distance_from_structure_file(
            file="../autodata/align/align_seed_sequences_with_structure/3rod.pdb",
            residue_dictionary={11:"Y"})

        list1=acs_dictionary[11]
        self.assertEqual(9,list1[0][0])
        self.assertEqual('ASP', list1[0][1])
        self.assertEqual(int(4.878653-list1[0][2]),0)


def main():
    unittest.main()
    print(0)

if __name__ == "__main__":
    main()