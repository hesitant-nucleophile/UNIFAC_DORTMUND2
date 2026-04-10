#put all the files in a folder, then run:

from folder_name import Unifac_Dortmund

smiles = ["CCO", "CC(=O)C"]
mol_frac_lst=[0.5, 0.5]
T=298
system=Unifac_Dortmund(smiles_lst,mol_frac_lst,T)

# gives activity coefficient of all components in the mixture
system.gamma_total()
# gives activity coefficient of component i in the mixture, according to the  order in which the components were passed.
system.gamma_singular(i)
