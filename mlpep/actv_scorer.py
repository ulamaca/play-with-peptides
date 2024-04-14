import torch
import pandas as pd
from mlpep.models import Classifier
from mlpep.datasetbioactivity import Dataset
from mlpep.datasetbioactivity import collate_fn_no_activity as collate_fn

model_path = '/home/ulamaca/projects/12_mlpeptide/MLpeptide/data/AIpep-clean/models/RNN-classifier/em100_hi400_la2_ep48'

def load_model(filename = model_path):        
    model = Classifier.load_from_file(filename)    

    return model

def predict_wo_label(data_loader, model):    
    out_list = []

    with torch.no_grad():
        for i_batch, sample_batched in enumerate(data_loader):                
            
            seq_batched = sample_batched[0].to(model.device, non_blocking=True)
            seq_lengths = sample_batched[1].to(model.device, non_blocking=True)            
            out_list += torch.exp(model.evaluate(seq_batched, seq_lengths))[: ,1].to("cpu", non_blocking=True)
        
        out_list = torch.stack(out_list)
    return out_list.cpu().numpy()

def get_pred_result(seqs: list, clf):
    '''
        prepare sequence into the format for peptide actv clf        
    '''
    df = pd.DataFrame([{'Sequence': seq} for seq in seqs])    
    gen_dataset = Dataset(df, clf.vocabulary, with_activity=False)    
    gen_dataloader = torch.utils.data.DataLoader(gen_dataset, batch_size=1, shuffle=False, collate_fn = collate_fn, drop_last=True, pin_memory=True, num_workers=4)
    y_score = predict_wo_label(gen_dataloader, clf)
    df['p_actv'] = y_score
    
    return df


if __name__ == "__main__":
    ## CLS model
    clf = load_model()