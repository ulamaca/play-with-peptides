
# IUPAC, standard AA chars (#=20, w/o B,J,O,U,X,Z)
normal_AA = ["A","C","D","E","F","G","H","I","L","M","N","P","K","Q","R","S","T","V","W","Y"]


# cite from https://github.com/reymond-group/PDGA-MAP4_AP/blob/master/pdga/sequence.py

# list of possible aminoacids (+ non-standard, + cycles etc. )
AA = ['R', 'H', 'K', 'E', 'S', 'T', 'N', 'Q', 'G', 'P', 'A', 'V', 'I', 'L', 'F', 'Y', 'W', 'C', 'D', 'M', 'Z', 'O',
        '!', '?', '=', '%', '$', '@', '#']

# map AA to SMILES
AA_SMILES = {'A': '[N:2][C@@H](C)[C:1](O)=O', 'R': '[N:2][C@@H](CCCNC(N)=N)[C:1](O)=O',
                'N': '[N:2][C@@H](CC(N)=O)[C:1](O)=O', 'D': '[N:2][C@@H](CC(O)=O)[C:1](O)=O',
                'C': '[N:2][C@@H](CS)[C:1](O)=O', 'Q': '[N:2][C@@H](CCC(N)=O)[C:1](O)=O',
                'E': '[N:2][C@@H](CCC(O)=O)[C:1](O)=O', 'G': '[N:2]C[C:1](O)=O',
                'H': '[N:2][C@@H](CC1=CNC=N1)[C:1](O)=O', 'I': '[N:2][C@@H]([C@@H](C)CC)[C:1](O)=O',
                'K': '[N:2][C@@H](CCCCN)[C:1](O)=O', 'L': '[N:2][C@@H](CC(C)C)[C:1](O)=O',
                'M': '[N:2][C@@H](CCSC)[C:1](O)=O', 'F': '[N:2][C@@H](CC1=CC=CC=C1)[C:1](O)=O',
                'P': 'C1CC[N:2][C@@H]1[C:1](O)=O', 'S': '[N:2][C@@H](CO)[C:1](O)=O',
                'T': '[N:2][C@@H]([C@H](O)C)[C:1](O)=O', 'W': '[N:2][C@@H](CC1=CNC2=CC=CC=C12)[C:1](O)=O',
                'Y': '[N:2][C@@H](CC1=CC=C(C=C1)O)[C:1](O)=O', 'V': '[N:2][C@@H](C(C)C)[C:1](O)=O',
                'Ä': '[N:2][C@@H](C[S:1])[C:1](O)=O', 'Ö': '[N:2][C@@H](C[S:2])[C:1](O)=O',
                'Ü': '[N:2][C@@H](C[S:3])[C:1](O)=O',
                'Z': 'C1C(O)C[N:2][C@@H]1[C:1](O)=O',
                'O': '[N:2][C@@H](CCCN)[C:1](O)=O',
                'a': '[N:2][C@H](C)[C:1](O)=O', 'r': '[N:2][C@H](CCCNC(N)=N)[C:1](O)=O',
                'n': '[N:2][C@H](CC(N)=O)[C:1](O)=O', 'd': '[N:2][C@H](CC(O)=O)[C:1](O)=O',
                'c': '[N:2][C@H](CS)[C:1](O)=O', 'q': '[N:2][C@H](CCC(N)=O)[C:1](O)=O',
                'e': '[N:2][C@H](CCC(O)=O)[C:1](O)=O', 'g': '[N:2]C[C:1](O)=O',
                'h': '[N:2][C@H](CC1=CNC=N1)[C:1](O)=O', 'i': '[N:2][C@H]([C@@H](C)CC)[C:1](O)=O',
                'k': '[N:2][C@H](CCCCN)[C:1](O)=O', 'l': '[N:2][C@H](CC(C)C)[C:1](O)=O',
                'm': '[N:2][C@H](CCSC)[C:1](O)=O', 'f': '[N:2][C@H](CC1=CC=CC=C1)[C:1](O)=O',
                'p': 'C1CC[N:2][C@H]1[C:1](O)=O', 's': '[N:2][C@H](CO)[C:1](O)=O',
                't': '[N:2][C@H]([C@H](O)C)[C:1](O)=O', 'w': '[N:2][C@H](CC1=CNC2=CC=CC=C12)[C:1](O)=O',
                'y': '[N:2][C@H](CC1=CC=C(C=C1)O)[C:1](O)=O', 'v': '[N:2][C@H](C(C)C)[C:1](O)=O',
                'ä': '[N:2][C@H](C[S:1])[C:1](O)=O', 'ö': '[N:2][C@H](C[S:2])[C:1](O)=O',
                'ü': '[N:2][C@H](C[S:3])[C:1](O)=O',
                '!': '[N:2]CC[C:1](O)=O', '?': '[N:2]CCC[C:1](O)=O',
                '=': '[N:2]CCCC[C:1](O)=O', '%': '[N:2]CCCCC[C:1](O)=O',
                '$': '[N:2]CCCCCC[C:1](O)=O', '@': '[N:2]CCCCCCC[C:1](O)=O',
                '#': '[N:2]CC[C:1](O)=O'}

# list of possible C-terminals
CT = ['+']
# list of possible N-capping
NT = ['&']

interprete_dict = {'Arg': 'R', 'His': 'H', 'Lys': 'K', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Thr': 'T', 'Asn': 'N',
                    'Gln': 'Q', 'Cys': 'C', 'Sec': 'U', 'Gly': 'G', 'Pro': 'P', 'Ala': 'A', 'Ile': 'I', 'Leu': 'L',
                    'Met': 'M', 'Phe': 'F', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Dap': '1', 'Dab': '2',
                    'BOrn': '3', 'BLys': '4', 'Hyp': 'Z', 'Orn': 'O', 'bAla': '!', 'Gaba': '?', 'dDap': '5',
                    'dDab': '6',
                    'dBOrn': '7', 'dBLys': '8', 'dArg': 'r', 'dHis': 'h', 'dLys': 'k', 'dAsp': 'd', 'dGlu': 'e',
                    'dSer': 's',
                    'dThr': 't', 'dAsn': 'n', 'dGln': 'q', 'dCys': 'c', 'dSec': 'u', 'dGly': 'g', 'dPro': 'p',
                    'dAla': 'a',
                    'dIle': 'i', 'dLeu': 'l', 'dMet': 'm', 'dPhe': 'f', 'dTrp': 'w', 'dTyr': 'y', 'dVal': 'v',
                    'dHyp': 'z', 'dOrn': 'o', 'a5a': '=', 'a6a': '%', 'a7a': '$', 'a8a': '@', 'a9a': '#',
                    'Cys1': 'Ä', 'Cys2': 'Ö', 'Cys3': 'Ü', 'dCys1': 'ä', 'dCys2': 'ö', 'dCys3': 'ü',
                    'Ac': '&', 'NH2': '+', 'met': '-', 'cy': 'X'}

# list of possible branching units (1=Dap, 2=Dab, 3=Orn, 4=Lys)
B = ['1', '2', '3', '4']
B_SMILES = {'1': '[N:2][C@@H](C[N:2])[C:1](O)=O', '2': '[N:2][C@@H](CC[N:2])[C:1](O)=O',
            '3': '[N:2][C@@H](CCC[N:2])[C:1](O)=O', '4': '[N:2][C@@H](CCCC[N:2])[C:1](O)=O',
            '5': '[N:2][C@H](C[N:2])[C:1](O)=O', '6': '[N:2][C@H](CC[N:2])[C:1](O)=O',
            '7': '[N:2][C@H](CCC[N:2])[C:1](O)=O', '8': '[N:2][C@H](CCCC[N:2])[C:1](O)=O'}