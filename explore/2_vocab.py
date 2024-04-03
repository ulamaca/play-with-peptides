import pandas as pd
import matplotlib.pyplot as plt
import json

if __name__ == "__main__":
    df_1 = pd.read_csv('./data/mlpeptide/DAASP_RNN_dataset.csv', sep=',', index_col='ID')
    df_2 = pd.read_csv('./data/mlpeptide/DBAASP_nointrabond.csv', sep=';', index_col='ID')

    df = pd.concat([df_1, df_2], axis=0)
    df.dropna(subset=['Sequence'], inplace=True)
    df = df.reset_index()

    char_counts = {}

    # Iterate over each sequence in the DataFrame
    for sequence in df['Sequence']:
        for char in sequence:
            # Increment the character count in the dictionary
            if char in char_counts:
                char_counts[char] += 1
            else:
                char_counts[char] = 1

    # char_counts now contains the counts of each unique character across all sequences
    print("the vocab and the statistics of Peptide Sequences:")
    print("#1 Counts:")
    char_counts = dict(sorted(char_counts.items(), key=lambda x: x[1], reverse=True))
    print(json.dumps(char_counts, indent=4))
    print("#1 Percentage:")
    total_counts = sum(list(char_counts.values()))
    char_percentage = {k:v/total_counts for k,v in char_counts.items()}
    print(json.dumps(char_percentage, indent=4))