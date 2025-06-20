"""test the GED workflow itself using the ChemGED package"""

import numpy as np
import pytest

from chemged.chemged import ApproximateChemicalGED
from chemged.cost import UniformElementCostMatrix
from chemged.plotting import plot_assignment


class TestGEDIntegration:
    """Integration tests for the ChemGED package."""

    @pytest.fixture
    def ged_calculator(self):
        """Create an ApproximateChemicalGED calculator with default cost matrix."""
        cost_matrix = UniformElementCostMatrix()
        return ApproximateChemicalGED(cost_matrix)

    def test_full_workflow(self, ged_calculator):
        """Test the full workflow from SMILES to GED calculation and visualization."""
        # Define some SMILES strings for testing
        smiles1 = "CC"  # Ethane
        smiles2 = "CO"  # Methanol

        # Calculate GED with assignment
        dist, assignment = ged_calculator.compute_ged(smiles1, smiles2, return_assignment=True)

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

        # Check that assignment is a tuple of two numpy arrays
        assert isinstance(assignment, tuple)
        assert len(assignment) == 2
        assert isinstance(assignment[0], np.ndarray)
        assert isinstance(assignment[1], np.ndarray)

        # Create a visualization (with show=False to avoid opening a window during tests)
        fig = plot_assignment(smiles1, smiles2, assignment, show=False)

        # Check that the result is a matplotlib Figure
        assert fig is not None

    def test_multiple_molecules(self, ged_calculator):
        """Test calculating GED between multiple molecules."""
        # Define some SMILES strings for testing
        smiles_list = ["CC", "CO", "CCO", "c1ccccc1"]  # Ethane, Methanol, Ethanol, Benzene

        # Calculate pairwise distances
        pdist = ged_calculator.pdist(smiles_list)

        # Check that the result has the correct shape
        n = len(smiles_list)
        expected_length = n * (n - 1) // 2
        assert pdist.shape == (expected_length,)

        # Calculate cross distances
        cdist = ged_calculator.cdist(smiles_list[:2], smiles_list[2:])

        # Check that the result has the correct shape
        n1 = 2
        n2 = 2
        expected_length = n1 * n2
        assert cdist.shape == (expected_length,)

    def test_custom_cost_matrix(self):
        """Test using a custom cost matrix."""
        # Create a custom cost matrix
        custom_matrix = UniformElementCostMatrix(
            node_sub_cost=2.0,
            node_ins_cost=3.0,
            node_del_cost=4.0,
            edge_sub_cost=5.0,
            edge_ins_cost=6.0,
            edge_del_cost=7.0,
        )

        # Create a GED calculator with the custom cost matrix
        ged_calculator = ApproximateChemicalGED(custom_matrix)

        # Define some SMILES strings for testing
        smiles1 = "CC"  # Ethane
        smiles2 = "CO"  # Methanol

        # Calculate GED
        dist = ged_calculator.compute_ged(smiles1, smiles2)

        # Check that the distance is a non-negative finite number
        assert dist >= 0
        assert np.isfinite(dist)

        # The distance should be different from the default cost matrix
        default_matrix = UniformElementCostMatrix()
        default_calculator = ApproximateChemicalGED(default_matrix)
        default_dist = default_calculator.compute_ged(smiles1, smiles2)

        # The custom costs are higher, so the distance should be higher
        assert dist > default_dist

    def test_triangle_inequality_approximate(self, ged_calculator):
        """
        Test that the GED approximately satisfies the triangle inequality.

        Note: The GED is an approximation, so it might not strictly satisfy the
        triangle inequality, but it should be close.
        """
        # Define some SMILES strings for testing
        smiles1 = "CC"  # Ethane
        smiles2 = "CO"  # Methanol
        smiles3 = "CCO"  # Ethanol

        # Calculate pairwise distances
        dist12 = ged_calculator.compute_ged(smiles1, smiles2)
        dist23 = ged_calculator.compute_ged(smiles2, smiles3)
        dist13 = ged_calculator.compute_ged(smiles1, smiles3)

        # Check that the triangle inequality approximately holds
        # d(a,c) <= d(a,b) + d(b,c)
        # Allow for some numerical error
        assert dist13 <= dist12 + dist23 + 1e-10
