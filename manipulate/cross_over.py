import random

def naive_cross_over(seq_1, seq_2, short_tolerance=3):
    '''
        seq_1: str
        seq_2: str
        short_tolerance: int
            how much #char shorter could one tolerate after cross over

        a naive way to cross-over two seq
    
    '''
    l_1 = len(seq_1)
    l_2 = len(seq_2)

    L = min(l_1, l_2)
    
    segment_pos = random.randint(1,L-2) # [<start>,x,y,z...,<end>], the seg should not be start nor end
    final_len = random.randint(L-short_tolerance, l_1+l_2)
    

    if random.uniform(0,1) < 0.5:
        seq = seq_1[:segment_pos] + seq_2[segment_pos:]
    else:
        seq = seq_2[:segment_pos] + seq_1[segment_pos:]

    return seq[:final_len]
