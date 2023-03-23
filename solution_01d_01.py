from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

theobromine = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)NC(=O)N2C') # insert SMILES
xanthine = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=O)N2') #insert SMILES
fig1 = Draw.MolToImage(theobromine, size=(500, 500))
fig1.save('./theo.png')
fig2 = Draw.MolToImage(xanthine, size=(500, 500))
fig2.save('./xant.png')

mols = [theobromine, xanthine] #create a list containing the 3 mol objects we have created
names = [ 'theobromine', 'xanthine'] #create a list containing the names of the 3 molecules

#Now we create the GridImage
grid = Draw.MolsToGridImage(mols, legends=names) #pass the 'mols' list here and create the image
