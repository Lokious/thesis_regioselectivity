import unittest

import pandas as pd
from Model_class import Model_class
class TestModel_class(unittest.TestCase):

    model = Model_class()
    def test0_check_file_exist(self):
        self.model.check_file_exist()
    def test1_keep_methyled_substrate(self):
       with self.assertRaises(AttributeError):
           self.model.keep_methyled_substrate(1)
    def test2_return_reactions(self):
        dataframe = pd.read_csv("../data/seq_smiles_all.csv",header=0,index_col=0)
        #the first reactant site should be O49:18
        self.assertEqual(dataframe.loc[0,"reactant_site"],"O49:18")
    def test3_return_reactions(self):
        dataframe = pd.read_csv("../data/seq_smiles_all.csv",header=0,index_col=0)
        #the first reactant site should be O49:18
        main_product = dataframe.loc[0,"main_sub"]
        main_substrate = dataframe.loc[0,"main_pro"]
        self.assertLess(main_substrate,main_product)



def main():
    unittest.main()
    print(0)

if __name__ == "__main__":
    main()
