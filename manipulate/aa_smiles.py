# cite from https://github.com/learningmatter-mit/peptimizer/blob/665fa51e89b98e8ef280a95c6f26a82ecdb57f18/dataset/data_cpp/cpp_smiles.json
aa2smiles_dict = {"2": "NC(CSC1=C(F)C(F)=C(C(F)=C1F)C1=C(F)C(F)=C(SCC(N)C(N)=O)C(F)=C1F)C(N)=O", 
                   "3": "CC(=O)CC1=CN(CCCCC(N)C(N)=O)N=N1", 
                   "A": "N[C@@H](C)C(O)=O", 
                   "B": "C(CN)C(=O)O", 
                   "X": "C(CCC(=O)O)CCN", 
                   "R": "N[C@@H](CCCNC(N)=N)C(O)=O", 
                   "N": "N[C@@H](CC(N)=O)C(O)=O", 
                   "D": "N[C@@H](CC(O)=O)C(O)=O", 
                   "C": "N[C@H](C(O)=O)CS", 
                   "E": "N[C@@H](CCC(O)=O)C(O)=O", 
                   "Q": "N[C@@H](CCC(N)=O)C(O)=O", 
                   "G": "NCC(O)=O", 
                   "H": "N[C@@H](CC1=CNC=N1)C(O)=O", 
                   "I": "N[C@@H]([C@@H](C)CC)C(O)=O", 
                   "L": "N[C@@H](CC(C)C)C(O)=O", 
                   "K": "N[C@@H](CCCCN)C(O)=O", 
                   "M": "N[C@@H](CCSC)C(O)=O", 
                   "F": "N[C@@H](CC1=CC=CC=C1)C(O)=O", 
                   "P": "O=C(O)[C@H]1NCCC1", 
                   "S": "N[C@@H](CO)C(O)=O", 
                   "T": "N[C@@H]([C@H](O)C)C(O)=O", 
                   "W": "N[C@@H](CC1=CNC2=C1C=CC=C2)C(O)=O", 
                   "Y": "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O", 
                   "V": "N[C@@H](C(C)C)C(O)=O", 
                   "@": "N[C@@H](CSC1=C(C(F)=C(C(F)=C1F)C2=C(C(F)=C(C(F)=C2F)SC[C@@H](C(O)=O)N)F)F)C(O)=O", 
                   "#": "N[C@H](C(O)=O)CSC1=CC(SC[C@@H](N)C(O)=O)=CC(SC[C@H](N)C(O)=O)=C1"}

# the second
#aa2simles_dict = {"A": "N[C@@]([H])(C)C(O)=O", "R": "N[C@@]([H])(CCCNC(NS(C1=C(C)C(CC(C)(C)O2)=C2C(C)=C1C)(=O)=O)=N)C(O)=O", "N": "N[C@@]([H])(CCC(NC(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3)=O)C(O)=O", "D": "N[C@@]([H])(CC(OC(C)(C)C)=O)C(O)=O", "C": "N[C@@]([H])(CSC(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3)C(O)=O", "E": "N[C@@]([H])(CCC(OC(C)(C)C)=O)C(O)=O", "Q": "N[C@@]([H])(CC(NC(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3)=O)C(O)=O", "G": "NCC(O)=O", "H": "N[C@@]([H])(CC1=CN=CN1C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C(O)=O", "I": "N[C@@]([H])([C@@H](C)CC)C(O)=O", "L": "N[C@@]([H])(CC(C)C)C(O)=O", "K": "N[C@@]([H])(CCCCNC(OC(C)(C)C)=O)C(O)=O", "M": "N[C@@]([H])(CCSC)C(O)=O", "F": "N[C@@]([H])(CC1=CC=CC=C1)C(O)=O", "P": "OC([C@]1([H])CCCN1)=O", "S": "N[C@@]([H])(COC(C)(C)C)C(O)=O", "T": "N[C@@]([H])([C@H](OC(C)(C)C)C)C(O)=O", "W": "N[C@@]([H])(CC1=CN(C(OC(C)(C)C)=O)C2=C1C=CC=C2)C(O)=O", "Y": "N[C@@]([H])(CC1=CC=C(OC(C)(C)C)C=C1)C(O)=O", "V": "N[C@@]([H])(C(C)C)C(O)=O"}