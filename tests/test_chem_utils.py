"""tests for chem_utils.py"""

import pytest
from networkx import Graph
from rdkit import Chem

from chemged.chem_utils import mol_to_nx, nx_to_mol, to_mol


class TestToMol:
    """Tests for the to_mol function in chem_utils.py."""

    def test_to_mol_with_mol_object(self):
        """Test that to_mol returns the same Mol object when given a Mol object."""
        mol = Chem.MolFromSmiles("CC")
        result = to_mol(mol, fail_on_error=True)
        assert result is mol

    def test_to_mol_with_valid_smiles(self):
        """Test that to_mol correctly converts a valid SMILES string to a Mol object."""
        smiles = "CC"
        result = to_mol(smiles, fail_on_error=True)
        assert isinstance(result, Chem.Mol)
        assert Chem.MolToSmiles(result) == smiles

    def test_to_mol_with_invalid_smiles_fail_on_error_true(self):
        """Test that to_mol raises ValueError when given an invalid SMILES string"""
        with pytest.raises(ValueError):
            to_mol("invalid_smiles", fail_on_error=True)

    def test_to_mol_with_invalid_smiles_fail_on_error_false(self):
        """Test that to_mol returns None when given an invalid SMILES string"""
        result = to_mol("invalid_smiles", fail_on_error=False)
        assert result is None

    def test_to_mol_with_invalid_type_fail_on_error_true(self):
        """Test that to_mol raises TypeError when given an invalid type"""
        with pytest.raises(TypeError):
            to_mol(123, fail_on_error=True)

    def test_to_mol_with_invalid_type_fail_on_error_false(self):
        """Test that to_mol returns None when given an invalid type"""
        result = to_mol(123, fail_on_error=False)
        assert result is None


class TestMolToNx:
    """Tests for the mol_to_nx function in chem_utils.py."""

    def test_mol_to_nx_simple_molecule(self):
        """Test that mol_to_nx correctly converts a simple molecule to a NetworkX graph."""
        mol = Chem.MolFromSmiles("CC")
        graph = mol_to_nx(mol)

        assert isinstance(graph, Graph)

        # Check that the graph has the correct number of nodes and edges
        assert len(graph.nodes) == 2
        assert len(graph.edges) == 1

        # Check that the nodes have the correct attributes
        for node in graph.nodes:
            assert "atomic_num" in graph.nodes[node]
            assert "hybridization" in graph.nodes[node]
            assert "is_aromatic" in graph.nodes[node]
            assert "is_in_ring" in graph.nodes[node]
            assert "hydrogen" in graph.nodes[node]
            assert "degree" in graph.nodes[node]

        # Check that the edges have the correct attributes
        for edge in graph.edges:
            assert "bond_type" in graph.get_edge_data(*edge)

    def test_mol_to_nx_complex_molecule(self):
        """Test that mol_to_nx correctly converts a more complex molecule to a NetworkX graph."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        graph = mol_to_nx(mol)

        # Check that the graph has the correct number of nodes and edges
        assert len(graph.nodes) == 6
        assert len(graph.edges) == 6

        # Check that all nodes are aromatic
        for node in graph.nodes:
            assert graph.nodes[node]["is_aromatic"] is True
            assert graph.nodes[node]["is_in_ring"] is True


class TestNxToMol:
    """Tests for the nx_to_mol function in chem_utils.py."""

    def test_nx_to_mol_simple_graph(self):
        """Test that nx_to_mol correctly converts a simple graph back to a Mol object."""
        # Create a molecule, convert to graph, then back to molecule
        original_mol = Chem.MolFromSmiles("CC")
        graph = mol_to_nx(original_mol)
        new_mol = nx_to_mol(graph)

        # Check that the new molecule has the same number of atoms and bonds
        assert new_mol.GetNumAtoms() == original_mol.GetNumAtoms()
        assert new_mol.GetNumBonds() == original_mol.GetNumBonds()

    def test_roundtrip_conversion(self):
        """Test that converting a molecule to a graph and back preserves the structure."""
        # Test with a more complex molecule
        smiles_list = ["c1ccccc1", "CC(=O)O", "C1=CC=C(C=C1)C(=O)O"]

        for smiles in smiles_list:
            original_mol = Chem.MolFromSmiles(smiles)
            graph = mol_to_nx(original_mol)
            new_mol = nx_to_mol(graph)

            # Check that the canonical SMILES of the original and new molecules match
            # Note: This might not always work due to information loss in the conversion
            # but should work for these simple cases
            assert Chem.MolToSmiles(new_mol) == Chem.MolToSmiles(original_mol)
