"""test for cost matrix for chemical graph edit distance"""

import numpy as np
import pytest
from rdkit import Chem

from chemged.chem_utils import mol_to_nx
from chemged.cost import UniformElementCostMatrix


class TestUniformElementCostMatrix:
    """Test suite for UniformElementCostMatrix."""

    @pytest.fixture
    def cost_matrix(self):
        """Create a UniformElementCostMatrix with default values."""
        return UniformElementCostMatrix()

    @pytest.fixture
    def custom_cost_matrix(self):
        """Create a UniformElementCostMatrix with custom values."""
        return UniformElementCostMatrix(
            node_sub_cost=2.0,
            node_ins_cost=3.0,
            node_del_cost=4.0,
            edge_sub_cost=5.0,
            edge_ins_cost=6.0,
            edge_del_cost=7.0,
        )

    @pytest.fixture
    def simple_graphs(self):
        """Create two simple molecular graphs for testing."""
        mol1 = Chem.MolFromSmiles("CC")  # Ethane
        mol2 = Chem.MolFromSmiles("CO")  # Methanol
        g1 = mol_to_nx(mol1)
        g2 = mol_to_nx(mol2)
        return g1, g2

    def test_initialization(self, cost_matrix, custom_cost_matrix):
        """Test that the cost matrix is initialized with the correct values."""
        # Default values
        assert cost_matrix.node_sub_cost == 1.0
        assert cost_matrix.node_ins_cost == 1.0
        assert cost_matrix.node_del_cost == 1.0
        assert cost_matrix.edge_sub_cost == 1.0
        assert cost_matrix.edge_ins_cost == 1.0
        assert cost_matrix.edge_del_cost == 1.0

        # Custom values
        assert custom_cost_matrix.node_sub_cost == 2.0
        assert custom_cost_matrix.node_ins_cost == 3.0
        assert custom_cost_matrix.node_del_cost == 4.0
        assert custom_cost_matrix.edge_sub_cost == 5.0
        assert custom_cost_matrix.edge_ins_cost == 6.0
        assert custom_cost_matrix.edge_del_cost == 7.0

    def test_get_node_substitution_costs(self, cost_matrix, simple_graphs):
        """Test that get_node_substitution_costs returns the correct cost matrix."""
        g1, g2 = simple_graphs
        cost_mat = cost_matrix.get_node_substitution_costs(g1, g2)

        # Check shape
        assert cost_mat.shape == (len(g1.nodes), len(g2.nodes))

        # Check values - should be 1.0 where atomic numbers are equal, 0.0 otherwise
        # Ethane has 2 carbon atoms, Methanol has 1 carbon and 1 oxygen
        # So the cost matrix should have a 1.0 at (0,0) and (1,0) for C-C substitutions
        # and 0.0 at (0,1) and (1,1) for C-O substitutions
        atomic_nums1 = [g1.nodes[n]["atomic_num"] for n in g1.nodes]
        atomic_nums2 = [g2.nodes[n]["atomic_num"] for n in g2.nodes]

        for i, a1 in enumerate(atomic_nums1):
            for j, a2 in enumerate(atomic_nums2):
                if a1 == a2:
                    assert cost_mat[i, j] == 0.0
                else:
                    assert cost_mat[i, j] == 1.0

    def test_get_node_insertion_costs(self, cost_matrix, simple_graphs):
        """Test that get_node_insertion_costs returns the correct cost vector."""
        g1, _ = simple_graphs
        costs = cost_matrix.get_node_insertion_costs(g1)

        # Check shape
        assert costs.shape == (len(g1.nodes),)

        # Check values - should be all 1.0
        assert np.all(costs == 1.0)

    def test_get_node_deletion_costs(self, cost_matrix, simple_graphs):
        """Test that get_node_deletion_costs returns the correct cost vector."""
        g1, _ = simple_graphs
        costs = cost_matrix.get_node_deletion_costs(g1)

        # Check shape
        assert costs.shape == (len(g1.nodes),)

        # Check values - should be all 1.0
        assert np.all(costs == 1.0)

    def test_get_edge_substitution_costs(self, cost_matrix, simple_graphs):
        """Test that get_edge_substitution_costs returns the correct cost matrix."""
        g1, g2 = simple_graphs

        # Get adjacency views for the first node in each graph
        n1 = g1[0]  # First node in ethane
        n2 = g2[0]  # First node in methanol

        cost_mat = cost_matrix.get_edge_substitution_costs(n1, n2)

        # Check shape
        assert cost_mat.shape == (len(n1), len(n2))

        # Check values - should be 1.0 where bond types are equal, 0.0 otherwise
        # Both ethane and methanol have single bonds, so the cost should be 1.0
        bond_types1 = [e["bond_type"] for e in n1.values()]
        bond_types2 = [e["bond_type"] for e in n2.values()]

        for i, b1 in enumerate(bond_types1):
            for j, b2 in enumerate(bond_types2):
                if b1 == b2:
                    assert cost_mat[i, j] == 0.0
                else:
                    assert cost_mat[i, j] == 1.0

    def test_get_edge_insertion_costs(self, cost_matrix, simple_graphs):
        """Test that get_edge_insertion_costs returns the correct cost vector."""
        g1, _ = simple_graphs
        n1 = g1[0]  # First node in ethane

        costs = cost_matrix.get_edge_insertion_costs(n1)

        # Check shape
        assert costs.shape == (len(n1),)

        # Check values - should be all 1.0
        assert np.all(costs == 1.0)

    def test_get_edge_deletion_costs(self, cost_matrix, simple_graphs):
        """Test that get_edge_deletion_costs returns the correct cost vector."""
        g1, _ = simple_graphs
        n1 = g1[0]  # First node in ethane

        costs = cost_matrix.get_edge_deletion_costs(n1)

        # Check shape
        assert costs.shape == (len(n1),)

        # Check values - should be all 1.0
        assert np.all(costs == 1.0)

    def test_custom_costs(self, custom_cost_matrix, simple_graphs):
        """Test that the custom cost matrix returns the correct values."""
        g1, g2 = simple_graphs

        # Test node costs
        node_sub_costs = custom_cost_matrix.get_node_substitution_costs(g1, g2)
        node_ins_costs = custom_cost_matrix.get_node_insertion_costs(g1)
        node_del_costs = custom_cost_matrix.get_node_deletion_costs(g2)

        # Check that substitution costs use the custom value
        atomic_nums1 = [g1.nodes[n]["atomic_num"] for n in g1.nodes]
        atomic_nums2 = [g2.nodes[n]["atomic_num"] for n in g2.nodes]

        for i, a1 in enumerate(atomic_nums1):
            for j, a2 in enumerate(atomic_nums2):
                if a1 != a2:
                    assert node_sub_costs[i, j] == 2.0

        # Check that insertion and deletion costs use the custom values
        assert np.all(node_ins_costs == 3.0)
        assert np.all(node_del_costs == 4.0)

        # Test edge costs
        n1 = g1[0]
        n2 = g2[0]

        edge_sub_costs = custom_cost_matrix.get_edge_substitution_costs(n1, n2)
        edge_ins_costs = custom_cost_matrix.get_edge_insertion_costs(n1)
        edge_del_costs = custom_cost_matrix.get_edge_deletion_costs(n2)

        # Check that substitution costs use the custom value
        bond_types1 = [e["bond_type"] for e in n1.values()]
        bond_types2 = [e["bond_type"] for e in n2.values()]

        for i, b1 in enumerate(bond_types1):
            for j, b2 in enumerate(bond_types2):
                if b1 != b2:
                    assert edge_sub_costs[i, j] == 5.0

        # Check that insertion and deletion costs use the custom values
        assert np.all(edge_ins_costs == 6.0)
        assert np.all(edge_del_costs == 7.0)
