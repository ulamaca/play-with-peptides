import ast # may need ast.literal_eval
import pandas as pd

hemo_path = '/home/ulamaca/projects/10_peptide/data/dbaasp/fine_tune_hemolysis.csv'

def replace_inf_to_float_inf(str_):
    return str_.replace('inf', 'float("inf")')

if __name__ == "__main__":
    df_hemo = pd.read_csv(hemo_path, index_col='ID')
    df_hemo['len'] = df_hemo['Sequence'].map(len)
    ast.literal_eval(df_hemo['hemolysis'].iloc[0])

    df_hemo['hemolysis'] = df_hemo['hemolysis'].map(lambda x: eval(replace_inf_to_float_inf(x)))

    # only those having labels are left:
    df_hemo = df_hemo.query('isNotHemolytic >= 0')

    print("#1 label distribution:")
    print(df_hemo.groupby('isNotHemolytic').size())

    breakpoint()