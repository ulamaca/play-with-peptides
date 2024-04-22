import pandas as pd

path_actv = '/home/ulamaca/projects/10_peptide/data/dbaasp/fine_tune_activity.csv'


if __name__ == "__main__":
    df_actv = pd.read_csv(path_actv, index_col='ID')
    df_actv = df_actv.iloc[:, 1:] # remove the 1st col
    df_actv['len'] = df_actv['Sequence'].map(len)
    
    print('#1 #-labels:')
    print(df_actv.groupby('activity').size())
    print('#2 len distribution')
    print(df_actv['len'].describe())
    print('2.1 len==25 is around 82-83% quantile')
    print(df_actv.quantile(0.82))
    print(df_actv.quantile(0.83))

    breakpoint()