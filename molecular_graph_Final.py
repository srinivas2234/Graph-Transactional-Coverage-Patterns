import numpy as np
from dgl import DGLGraph
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdDesc
from utils import one_of_k_encoding_unk, one_of_k_encoding
import torch
import warnings

edge_list = []
vertex_list = []

def get_atom_features(atom, stereo, features, explicit_H=False):
    """
    Method that computes atom level features from rdkit atom object
    :param atom:
    :param stereo:
    :param features:
    :param explicit_H:
    :return: the node features of an atom
    """
    possible_atoms = ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I', 'Si']
    atom_features = one_of_k_encoding_unk(atom.GetSymbol(), possible_atoms)
    atom_features += one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1])
    atom_features += one_of_k_encoding_unk(atom.GetNumRadicalElectrons(), [0, 1])
    atom_features += one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6])
    atom_features += one_of_k_encoding_unk(atom.GetFormalCharge(), [-1, 0, 1])
    atom_features += one_of_k_encoding_unk(atom.GetHybridization(), [
        Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D])
    atom_features += [int(i) for i in list("{0:06b}".format(features))]

    if not explicit_H:
        atom_features += one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4])

    try:
        atom_features += one_of_k_encoding_unk(stereo, ['R', 'S'])
        atom_features += [atom.HasProp('_ChiralityPossible')]
    except Exception as e:

        atom_features += [False, False
                          ] + [atom.HasProp('_ChiralityPossible')]

    return np.array(atom_features)


def get_bond_features(bond):
    """
    Method that computes bond level features from rdkit bond object
    :param bond: rdkit bond object
    :return: bond features, 1d numpy array
    """

    bond_type = bond.GetBondType()
    bond_feats = [
        bond_type == Chem.rdchem.BondType.SINGLE, bond_type == Chem.rdchem.BondType.DOUBLE,
        bond_type == Chem.rdchem.BondType.TRIPLE, bond_type == Chem.rdchem.BondType.AROMATIC,
    ]
    # bond_feats += one_of_k_encoding_unk(str(bond.GetStereo()), ["STEREONONE", "STEREOANY", "STEREOZ", "STEREOE"])

    return np.array(bond_feats)


def get_graph_from_smile(molecule_smile):
    """
    Method that constructs a molecular graph with nodes being the atoms
    and bonds being the edges.
    :param molecule_smile: SMILE sequence
    :return: DGL graph object, Node features and Edge features
    """

    G = DGLGraph()
    molecule = Chem.MolFromSmiles(molecule_smile)
    features = rdDesc.GetFeatureInvariants(molecule)

    stereo = Chem.FindMolChiralCenters(molecule)
    chiral_centers = [0] * molecule.GetNumAtoms()
    for i in stereo:
        chiral_centers[i[0]] = i[1]

    G.add_nodes(molecule.GetNumAtoms())
    node_features = []
    edge_features = []
    for i in range(molecule.GetNumAtoms()):
        possible_atoms = ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I', 'Si']
        atom_i = molecule.GetAtomWithIdx(i)
        vertex_list.append([i,possible_atoms.index(atom_i.GetSymbol())])
        atom_i_features = get_atom_features(atom_i, chiral_centers[i], features[i])
        node_features.append(atom_i_features)

        for j in range(molecule.GetNumAtoms()):
            bond_ij = molecule.GetBondBetweenAtoms(i, j)
            if bond_ij is not None:
                G.add_edge(i, j)
                
                bond_typex = 0
                if(bond_ij.GetBondType() ==  Chem.rdchem.BondType.SINGLE):
                    bond_typex = 0

                elif(bond_ij.GetBondType() ==  Chem.rdchem.BondType.DOUBLE):
                    bond_typex = 1

                elif(bond_ij.GetBondType() ==  Chem.rdchem.BondType.TRIPLE):
                    bond_typex = 2

                elif(bond_ij.GetBondType() ==  Chem.rdchem.BondType.AROMATIC):
                    bond_typex = 3

                if(edge_list.count([j, i, bond_typex]) == 0):
                	edge_list.append([i,j,bond_typex])
                bond_features_ij = get_bond_features(bond_ij)
                edge_features.append(bond_features_ij)

    G.ndata['x'] = torch.FloatTensor(node_features)
    G.edata['w'] = torch.FloatTensor(edge_features)
    return G

file1 = open("./Dataset/Smilesdata.txt", "r")
file2 = open("./Dataset/Smilesdata_to_GT.txt", "w")
countx = 0
while True:
 
    # Get next line from file
    line = file1.readline()
 
    # if line is empty
    # end of file is reached
    if not line:
        break
    
    vertex_list.clear()
    edge_list.clear()
    smiles = line.split(',')
    file2.write("t # " + str(countx) + "\n")
    get_graph_from_smile(smiles[0])
    
    for v in vertex_list:
        file2.write("v " + str(v[0]) + " " + str(v[1]) + "\n")

    for e in edge_list:
        file2.write("e " + str(e[0]) + " " + str(e[1]) + " " + str(e[2]) + "\n")

    countx += 1


file1.close()
file2.close()
