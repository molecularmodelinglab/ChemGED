"""Test for ApproximateChemicalGED class in chemged module."""

import numpy as np
import pytest
from rdkit import Chem

from chemged.chem_utils import mol_to_nx
from chemged.chemged import MAX_INT32, ApproximateChemicalGED
from chemged.cost import UniformElementCostMatrix


class TestApproximateChemicalGED:
    """Test suite for ApproximateChemicalGED class."""

    @pytest.fixture
    def ged_calculator(self):
        """Create an ApproximateChemicalGED calculator with default cost matrix."""
        cost_matrix = UniformElementCostMatrix()
        return ApproximateChemicalGED(cost_matrix)

    @pytest.fixture
    def molecules(self):
        """Create simple molecules for testing."""
        mol1 = Chem.MolFromSmiles(
            "CCc(c1)ccc2[n+]1ccc3c2[nH]c4c3cccc4CCc1c[n+]2ccc3c4ccccc4[nH]c3c2cc1"
        )
        mol2 = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        return mol1, mol2

    @pytest.fixture
    def graphs(self, molecules):
        """Convert simple molecules to NetworkX graphs."""
        mol1, mol2 = molecules
        g1 = mol_to_nx(mol1)
        g2 = mol_to_nx(mol2)
        return g1, g2

    def test_initialization(self, ged_calculator):
        """Test that the GED calculator is initialized with the correct cost matrix."""
        assert isinstance(ged_calculator.cost_matrix, UniformElementCostMatrix)

    def test_edge_cost_matrix(self, ged_calculator, graphs):
        """Test that _edge_cost_matrix returns a matrix of the correct shape and values."""
        g1, g2 = graphs

        # Get adjacency views for nodes in each graph
        n1 = g1[2]
        n2 = g2[4]

        cost_mat = ged_calculator._edge_cost_matrix(n1, n2)

        # Check shape
        assert cost_mat.shape == (len(n1) + len(n2), len(n1) + len(n2))

        # check that the insert section has the correct structure
        for i in range(len(n2)):
            assert cost_mat[len(n1) + i, i] < MAX_INT32
            cost_mat[len(n1) + i, i] = MAX_INT32
        assert np.all(cost_mat[len(n1) :, 0 : len(n2)] == MAX_INT32)

        # Diagonal of deletion section should have finite values
        for i in range(len(n1)):
            assert cost_mat[i, len(n2) + i] < MAX_INT32
            cost_mat[i, len(n2) + i] = MAX_INT32
        assert np.all(cost_mat[0 : len(n1), len(n2) :] == MAX_INT32)

        # Substitution section (top left)
        assert np.all((cost_mat[0 : len(n1), 0 : len(n2)]) != MAX_INT32)

    def test_edge_edit_cost(self, ged_calculator, graphs):
        """Test that _edge_edit_cost returns a valid cost."""
        g1, g2 = graphs

        # Get adjacency views for the first node in each graph
        n1 = g1[0]
        n2 = g2[0]

        cost = ged_calculator._edge_edit_cost(n1, n2)

        # Check that the cost is a non-negative finite number
        assert cost >= 0
        assert np.isfinite(cost)

    def test_get_full_cost_matrix(self, ged_calculator, graphs):
        """Test that _get_full_cost_matrix returns a matrix of the correct shape and values."""
        g1, g2 = graphs

        cost_mat = ged_calculator._get_full_cost_matrix(g1, g2)

        # Check shape
        assert cost_mat.shape == (len(g1) + len(g2), len(g1) + len(g2))

        # Check that the matrix has the correct structure
        # Insertion section (bottom left)
        for i in range(len(g2)):
            assert cost_mat[len(g1) + i, i] < MAX_INT32
            cost_mat[len(g1) + i, i] = MAX_INT32
        assert np.all(cost_mat[len(g1) :, 0 : len(g2)] == MAX_INT32)

        # Deletion section (top right)
        for i in range(len(g1)):
            assert cost_mat[i, len(g2) + i] < MAX_INT32
            cost_mat[i, len(g2) + i] = MAX_INT32
        assert np.all(cost_mat[0 : len(g1), len(g2) :] == MAX_INT32)

        # Substitution section (top left)
        assert np.all((cost_mat[0 : len(g1), 0 : len(g2)]) < MAX_INT32)

    def test_ged_without_assignment(self, ged_calculator, graphs):
        """Test that _ged returns a valid distance without assignment."""
        g1, g2 = graphs

        dist, assignment = ged_calculator._ged(g1, g2, return_assignment=False)

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

        # Check that assignment is None
        assert assignment is None

    def test_ged_with_assignment(self, ged_calculator, graphs):
        """Test that _ged returns a valid distance with assignment."""
        g1, g2 = graphs

        dist, assignment = ged_calculator._ged(g1, g2, return_assignment=True)

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

        # Check that assignment is a tuple of two numpy arrays
        assert isinstance(assignment, tuple)
        assert len(assignment) == 2
        assert isinstance(assignment[0], np.ndarray)
        assert isinstance(assignment[1], np.ndarray)

        # Check that the assignment arrays have the same length
        assert len(assignment[0]) == len(assignment[1])

    def test_compute_ged_without_assignment(self, ged_calculator, molecules):
        """Test that compute_ged returns a valid distance without assignment."""
        mol1, mol2 = molecules

        dist = ged_calculator.compute_ged(mol1, mol2, return_assignment=False)

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

    def test_compute_ged_with_assignment(self, ged_calculator, molecules):
        """Test that compute_ged returns a valid distance with assignment."""
        mol1, mol2 = molecules

        result = ged_calculator.compute_ged(mol1, mol2, return_assignment=True)

        # Check that the result is a tuple
        assert isinstance(result, tuple)
        assert len(result) == 2

        dist, assignment = result

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

        # Check that assignment is a tuple of two numpy arrays
        assert isinstance(assignment, tuple)
        assert len(assignment) == 2
        assert isinstance(assignment[0], np.ndarray)
        assert isinstance(assignment[1], np.ndarray)

        # Check that the assignment arrays have the same length
        assert len(assignment[0]) == len(assignment[1])

    def test_compute_ged_with_smiles(self, ged_calculator):
        """Test that compute_ged works with SMILES strings."""
        smiles1 = "CC"  # Ethane
        smiles2 = "CO"  # Methanol

        dist = ged_calculator.compute_ged(smiles1, smiles2)

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

    def test_pdist(self, ged_calculator):
        """Test that pdist returns a vector of the correct shape and values."""
        chemicals = ["CC", "CO", "CCO"]  # Ethane, Methanol, Ethanol

        dists = ged_calculator.pdist(chemicals)

        # Check shape - should be (n*(n-1)/2,) where n is the number of chemicals
        n = len(chemicals)
        expected_length = n * (n - 1) // 2
        assert dists.shape == (expected_length,)

        # Check that all distances are non-negative finite numbers
        assert np.all(dists >= 0)
        assert np.all(np.isfinite(dists))

    def test_cdist(self, ged_calculator):
        """Test that cdist returns a vector of the correct shape and values."""
        chemicals1 = ["CC", "CO"]  # Ethane, Methanol
        chemicals2 = ["CCO", "c1ccccc1"]  # Ethanol, Benzene

        dists = ged_calculator.cdist(chemicals1, chemicals2)

        # Check shape - should be (n*m,) where n is the number of chemicals in chemicals1
        # and m is the number of chemicals in chemicals2
        n = len(chemicals1)
        m = len(chemicals2)
        expected_length = n * m
        assert dists.shape == (expected_length,)

        # Check that all distances are non-negative finite numbers
        assert np.all(dists >= 0)
        assert np.all(np.isfinite(dists))

    def test_identity_distance(self, ged_calculator):
        """Test that the distance between a molecule and itself is 0."""
        mol = Chem.MolFromSmiles("CC")

        dist = ged_calculator.compute_ged(mol, mol, return_assignment=False)

        assert dist == 0.0

    def test_symmetry(self, ged_calculator, molecules):
        """Test that the distance is symmetric."""
        mol1, mol2 = molecules

        dist1 = ged_calculator.compute_ged(mol1, mol2, return_assignment=False)
        dist2 = ged_calculator.compute_ged(mol2, mol1, return_assignment=False)

        assert dist1 == dist2
