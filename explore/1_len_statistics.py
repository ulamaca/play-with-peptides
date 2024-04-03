import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    df_1 = pd.read_csv('./data/mlpeptide/DAASP_RNN_dataset.csv', sep=',', index_col='ID')
    df_2 = pd.read_csv('./data/mlpeptide/DBAASP_nointrabond.csv', sep=';', index_col='ID')

    df = pd.concat([df_1, df_2], axis=0)
    # focus on EDA of df_1
    
    df['seq_len'] = df['Sequence'].apply(lambda x: len(x) if type(x)==str else None)
    df.dropna(subset=['Sequence'], inplace=True)
    df = df.reset_index()
    print(df['seq_len'].describe())
    df.hist('seq_len', bins=100)
    plt.show()
    #fig, ax = plt.figure()
    breakpoint()