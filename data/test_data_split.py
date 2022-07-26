from sklearn.model_selection import GroupShuffleSplit
import pandas as pd
data = pd.read_excel("test.xlsx", header = 0,index_col=0)
splitter = GroupShuffleSplit(test_size=.20, n_splits=1, random_state=7)
split = splitter.split(data, groups=df['molecular_id'])
train_inds, test_inds = next(split)

train = df.iloc[train_inds]
test = df.iloc[test_inds]
print(test,train)
