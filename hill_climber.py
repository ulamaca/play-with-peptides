import pandas as pd
from manipulate.mutate import Genetic_Mutations
from mlpep.actv_scorer import load_model, get_pred_result

def legalize_seq_for_clf(seq: str):
    '''
        as the title suggests
        TODO: to see why this vocab are not valid for the given vocab
    '''
    seq = seq.replace(' ', '')
    seq = seq.replace('X', '')
    seq = seq.replace('B', '')

    return seq


if __name__ == "__main__":
    # 1 define mutator
    data_path = 'data/cpp/cpp_predictor_dataset.csv'
    mutator = Genetic_Mutations(data_path=data_path)
    
    # 2 define pred model
    pred_model = load_model()
    
    # 3 hill climber
    n_step = 50
    n_population = 20

    seed_seq = 'TKPRPGP' # peptide: Selank

    generations = []
    
    # initial population
    population = [legalize_seq_for_clf(mutator.mutate(seed_seq)) for _ in range(n_population)]    

    for step in range(n_step):                
        # get prediction    
        df = get_pred_result(population, pred_model)
        df['generation'] = step + 1
        
        # generation trace
        generations.append(df)
                
        # mutate & legalize
        df_best = df.sort_values('p_actv').tail(1)
        best_step = df_best['Sequence'].item()
        score = df_best['p_actv'].item()
        population = [legalize_seq_for_clf(mutator.mutate(best_step)) for _ in range(n_population-1)] + [best_step]                

        print(f"step={step}, score={score}, current seed peptide: {best_step}")    

        
    df = pd.concat(generations).reset_index(drop=True)

    breakpoint()