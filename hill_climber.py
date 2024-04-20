'''
    240416: addded crossover to increase diversity
    240420: support a new cross_over: cut_inverted_cross_over
'''

import pandas as pd
import random
from manipulate.mutate import Genetic_Mutations
from manipulate.cross_over import naive_cross_over, cut_inverted_cross_over, random_inverted_cross_over
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
    n_population = 20 # real poputation size = n_population * (1.3) + 1
    # mutate_rate = 0.5
    cross_over_rate = 0.8
    cross_over_type = 'c_icv'
    random.seed(42)

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

        # TODO: to sample survivors (seq being intact) according to its current score

        # cross over
        mates = [random.sample(population, 2) for _ in range(int(n_population*cross_over_rate))]
        
        if cross_over_type == 'r_icv':
            offsprings = [random_inverted_cross_over(parent[0], parent[1]) for parent in mates]
        elif cross_over_type == 'c_icv':
            offsprings = [cut_inverted_cross_over(parent[0], parent[1]) for parent in mates]
        elif cross_over_type == 'naive':
            offsprings = [naive_cross_over(parent[0], parent[1]) for parent in mates]
        
        offsprings = [legalize_seq_for_clf(c) for c in offsprings]
        population.extend(offsprings)    
        population = list(set(population) ) # make them unique()

        print(f"step={step}, score={score}, current population size= {len(population)}, current seed peptide: {best_step}")    

        
    df = pd.concat(generations).reset_index(drop=True)

    breakpoint()