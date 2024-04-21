import pandas as pd

# df_2 is the mic dataset used for MIC classification at HyDRAMP paper

if __name__ == "__main__":
    mic_data_path = '/home/ulamaca/projects/10_peptide/data/hydramp/data/data/mic_data.csv'

    
    df_1 = pd.read_csv(mic_data_path, index_col=0)
    df_1['length'] = df_1['sequence'].apply(lambda x: len(x))
    
    df_2 = df_1.query('length <= 25')
    df_2['activity'] = df_2['value'].apply( lambda x: 'positive' if x <= 1.5 else 'negative')

    # 
    print("#1 the 'value' column distribution:")
    print(df_2.value.describe())

    print("#2 peptide sequence length distribution:")
    print(df_2.length.describe())

    print("#3 num. of positive vs negative data points")
    print(df_2.groupby('activity').size())

    breakpoint()