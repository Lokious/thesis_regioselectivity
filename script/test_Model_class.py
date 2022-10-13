import unittest

import pandas as pd
from Model_class import Model_class
class TestModel_class(unittest.TestCase):

    model = Model_class()
    def test0_check_file_exist(self):
        self.model.check_file_exist()

    def test1_keep_methyled_substrate(self):
       with self.assertRaises(RuntimeError):
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

    def test4_runPCA(self):
        """
        This function is to test if runPCA runs as expected

        """
        dataframe=pd.DataFrame({"A":[1,0,1],"B":[0,0,0],"methyl_type":["O","C","O"]},index=[0,1,2])

        y=pd.DataFrame(columns=["Y"],data=[1,0,1])
        pca_df=self.model.run_PCA(dataframe,y_label=y["Y"])
        self.assertEqual(0,pca_df.PC2[0])

    # def test5_create_similarity_matrix(self):
    #     self.model.create_similarity_matrix()

    def test6_use_less_similar_data_for_test(self):
        """
        This test take around 5 to 10 min
        :return:
        """
        input = pd.read_csv(
            "../autodata/input_data/input128fg_bond3_S_k_mer.csv", header=0,
            index_col=0)

        X_train, X_test, Y_train, Y_test=self.model.use_less_similar_data_for_test(input,num_bit=128)
        print(X_test)
        print(Y_test)
def main():
    unittest.main()
    print(0)

if __name__ == "__main__":
    main()
