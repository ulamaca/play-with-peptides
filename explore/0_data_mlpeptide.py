import pandas as pd

if __name__ == "__main__":
    df_1 = pd.read_csv('./data/mlpeptide/DAASP_RNN_dataset.csv', sep=',', index_col='ID')
    df_2 = pd.read_csv('./data/mlpeptide/DBAASP_nointrabond.csv', sep=';', index_col='ID')

    # focus on EDA of df_1
    df_1['seq_len'] = df_1['Sequence'].apply(len)
    print("Explorative Data Analysis:")
    print("###df_1 (RNN_dataset)")
    print("columns of df_1 are: ", df_1.columns, '\n')
    print(df_1.groupby('Set').size(), '\n')
    print(df_1.groupby('activity').size(), '\n')
    print(df_1['seq_len'].describe(), '\n')
    print("###")

    breakpoint()