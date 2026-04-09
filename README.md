put all the files in a folder, then run:

from folder_name import Unifac_Dortmund


smiles = ["CCO", "CC(=O)C"]
mol_lst=[0.5, 0.5]
T=298
system=Unifac_Dortmund(smiles_lst,mol_lst,T)
system.gamma_total()
