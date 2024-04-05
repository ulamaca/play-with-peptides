import pandas as pd


PATH_RNN = './data/mlpeptide/DAASP_RNN_dataset.csv'
PATH_2 = './data/mlpeptide/DBAASP_nointrabond.csv'


def get_rnn_dataset():
    return pd.read_csv(PATH_RNN, sep=',', index_col='ID')

def get_nib_dataset():
    return pd.read_csv(PATH_2, sep=';', index_col='ID')

def get_joint_dataset():
    df_1 = get_rnn_dataset()
    df_2 = get_nib_dataset()
    df = pd.concat([df_1, df_2], axis=0)
    df.dropna(subset=['Sequence'], inplace=True)
    df = df.reset_index()

    return df


