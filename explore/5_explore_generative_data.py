'''
    data for pretraining language models 
    (cited from Nature Comm. 2023 https://www.nature.com/articles/s41467-023-36994-z)
'''
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv('/home/ulamaca/projects/10_peptide/data/hydramp/data/data/Uniprot_0_25_train.csv')
    df['len'] = df['Sequence'].map(len)

    print("#0 columns")
    print(df.columns)
    print(df.head(5))
    print("#1 length statisitics")
    print(df.len.describe())

    breakpoint()
