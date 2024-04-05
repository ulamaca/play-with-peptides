import os
import sys
dir_ =  os.path.dirname(__file__) # the dir of the file
parent_dir = os.path.dirname(dir_)
sys.path.append(parent_dir)

from utils.data_dbaasp import get_joint_dataset, get_rnn_dataset
from utils.phys_chem_propterties import calculate_physchem_prop


if __name__ == "__main__":
    df = get_rnn_dataset()
    seqs = df.Sequence.to_list()

    # preprop
    # TODO: to know what does un-captical char mean in a peptide seq
    seqs = [s.upper() for s in seqs ] # the module can't hanlde non-capital characters
    pc_vals = calculate_physchem_prop(seqs)
    # ad pc_vals back into the df
    for k, vs in pc_vals.items():
        df[k] = vs    