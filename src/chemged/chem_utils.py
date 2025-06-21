"""Utilities for converting between RDKit Mol objects and NetworkX Graphs."""

from typing import Literal, Union, overload

import networkx as nx
from rdkit import Chem
from rdkit.rdBase import BlockLogs


Molable = Union[str, Chem.Mol]


@overload
def to_mol(smi: Molable, fail_on_error: Literal[False]) -> Chem.Mol: ...


@overload
def to_mol(smi: Molable, fail_on_error: Literal[True]) -> Chem.Mol | None: ...


def to_mol(smi: Molable, fail_on_error: bool = True) -> Chem.Mol | None:
    """
    Given an object, attempt to convert it to a Chem.Mol object

    Notes
    -----
    Can only covert valid SMILES str and rdkit.Mol objects.
    If a rdkit.Mol object is passed, the same object will be returned.

    Parameters
    ----------
    smi: Any
        object to convert to a Chem.Mol
    fail_on_error: bool
        whether to raise an exception when converting fails.
        if True, will return None when a conversion fails

    Returns
    -------
    Chem.Mol | None
        Will return a Chem.Mol object if conversion is successful,
        if fail_on_error is True and conversion fails, will return None.

    Raises
    ------
    ValueError
        if the SMILES cannot be parsed by rdkit
    TypeError
        if the passed object is not a type that can be converted to a Chem.Mol
    """
    _lock = BlockLogs()  # this turns off the rdkit logger
    if isinstance(smi, Chem.Mol):  # return if already a Chem.Mol
        return smi
    elif isinstance(smi, str):
        _mol = Chem.MolFromSmiles(smi)
        if _mol is None:
            if fail_on_error:
                raise ValueError(f"SMILES {smi} cannot be parsed by RDKit")
        return _mol
    else:
        if fail_on_error:
            raise TypeError(f"cannot convert type {type(smi)} to type rdkit.Mol")
        else:
            return None


def mol_to_nx(mol: Chem.Mol) -> nx.Graph:
    """
    Coverts an RDKit Mol object to a NetworkX Graph.

    Will store the following attributes for each node (atom):
    - atomic_num: Atomic number of the atom
    - hybridization: Hybridization state of the atom
    - is_aromatic: Whether the atom is aromatic
    - is_in_ring: Whether the atom is part of a ring
    - hydrogen: Total number of hydrogen atoms attached to the atom
    - degree: Degree of the atom (number of bonds)

    Will save only the bond type for each edge (bond)

    See the RDKit documentation for more details on the attributes,
    all of them are functionally Enums and can be treated as integers.

    Parameters
    ----------
    mol: Chem.Mol

    Returns
    -------
    nx.Graph
    """
    graph = nx.Graph()
    max_idx = 0
    for atom in mol.GetAtoms():
        graph.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            hybridization=atom.GetHybridization(),
            is_aromatic=atom.GetIsAromatic(),
            is_in_ring=atom.IsInRing(),
            hydrogen=atom.GetTotalNumHs(),
            degree=atom.GetDegree(),
            formal_charge=atom.GetFormalCharge(),
        )
        if atom.GetIdx() > max_idx:
            max_idx = atom.GetIdx()
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType())
    return graph


def nx_to_mol(graph: nx.Graph) -> Chem.Mol:
    """
    Give a NetworkX Graph, convert it to a Chem.Mol object.

    NOTE: This is only stable for converting graphs that were created from
    the `mol_to_nx` function in this package.
    It is also lossy, so some information may be lost in the conversion for
    Mol -> Graph -> Mol.

    Parameters
    ----------
    graph: nx.Graph
        the NetworkX Graph to convert to a Chem.Mol object

    Returns
    -------
    Chem.Mol
    """
    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(graph, "atomic_num")
    node_is_aromatics = nx.get_node_attributes(graph, "is_aromatic")
    node_hybridization = nx.get_node_attributes(graph, "hybridization")
    node_to_idx = {}
    for node in graph.nodes():
        a = Chem.Atom(atomic_nums[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridization[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    bond_types = nx.get_edge_attributes(graph, "bond_type")
    for edge in graph.edges():
        first, second = edge
        idx_first = node_to_idx[first]
        idx_second = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(idx_first, idx_second, bond_type)

    Chem.SanitizeMol(mol)
    return mol
