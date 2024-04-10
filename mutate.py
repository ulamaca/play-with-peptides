from manipulate.mutate import Genetic_Mutations


if __name__ == "__main__":
    sequence = 'TKPRPGP' # peptide: Selank
    data_path = 'data/cpp/cpp_predictor_dataset.csv'

    mutator = Genetic_Mutations(data_path=data_path)


    breakpoint()