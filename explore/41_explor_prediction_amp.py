import pandas as pd


def add_len_column(df: pd.DataFrame):
    assert 'Sequence' in df.columns
    df['length'] = df['Sequence'].apply(lambda x: len(x))

    return df

if __name__ == "__main__":
    amp_positive_data_path = '/home/ulamaca/projects/10_peptide/data/hydramp/data/data/unlabelled_positive.csv'
    amp_negative_data_path = '/home/ulamaca/projects/10_peptide/data/hydramp/data/data/unlabelled_negative.csv'

    df_p = pd.read_csv(amp_positive_data_path)
    df_p = add_len_column(df_p)
    df_p_1 = df_p.query('length <= 25')

    
    df_n = pd.read_csv(amp_negative_data_path)
    df_n = add_len_column(df_n)
    df_n_1 = df_n.query('length <= 25') # >> resulting empty df
    # in HyDramp (2023) paper, requires at least two steps to preprocess the negative dataset
        # 1 CD-HIT clustering
        # 2 croping the (longer) sequences to l <= 25
    
    print("#1 sequence length distribution of negative (raw) dataset", df_n.length.describe())

    breakpoint()
    