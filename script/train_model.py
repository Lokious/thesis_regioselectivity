import dill

from script import parse_data
from script.Model_class import Model_class
import pandas as pd
import numpy as np
import copy
import time
from datetime import date

mo_del = Model_class()
today = date.today()
# dd/mm/YY
d1 = today.strftime("%d_%m_%Y")
with open("data/input_data/input128fg_dpna_bi_{}".format(d1), 'rb') as file1:
    input_dataframe = dill.load(file1)

X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
    input_dataframe)

# mo_del.three_D_pca(X_train, y_train, "128fg")
mo_del.run_PCA(X_train, y_train, "128fg")
X_train = X_train.drop(columns=["methyl_type"])
X_test = X_test.drop(columns=["methyl_type"])
y_train = y_train.drop(columns=["methyl_type"])
y_test = y_test.drop(columns=["methyl_type"])
model = mo_del.RF_model(X_train, X_test, y_train, y_test,
                        "_input128fg_bi_type_{}".format(d1))