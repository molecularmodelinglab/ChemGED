"""test plotting functions for ApproximateChemicalGED"""

import os
import tempfile

import matplotlib.pyplot as plt
import pytest
from rdkit import Chem

from chemged.chemged import ApproximateChemicalGED
from chemged.cost import UniformElementCostMatrix
from chemged.plotting import plot_assignment


class TestPlotting:
    """test plotting utility"""

    @pytest.fixture
    def ged_calculator(self):
        """Create an ApproximateChemicalGED calculator with default cost matrix."""
        cost_matrix = UniformElementCostMatrix()
        return ApproximateChemicalGED(cost_matrix)

    @pytest.fixture
    def simple_molecules(self):
        """Create simple molecules for testing."""
        mol1 = Chem.MolFromSmiles("CC")  # Ethane
        mol2 = Chem.MolFromSmiles("CO")  # Methanol
        return mol1, mol2

    @pytest.fixture
    def assignment(self, ged_calculator, simple_molecules):
        """Get the assignment between two simple molecules."""
        mol1, mol2 = simple_molecules
        _, assignment = ged_calculator.compute_ged(mol1, mol2, return_assignment=True)
        return assignment

    def test_plot_assignment(self, simple_molecules, assignment):
        """Test that plot_assignment returns a matplotlib Figure."""
        mol1, mol2 = simple_molecules

        # Test with show=False to avoid opening a window during tests
        fig = plot_assignment(mol1, mol2, assignment, show=False)

        # Check that the result is a matplotlib Figure
        assert isinstance(fig, plt.Figure)

    def test_plot_assignment_with_smiles(self, assignment):
        """Test that plot_assignment works with SMILES strings."""
        smiles1 = "CC"  # Ethane
        smiles2 = "CO"  # Methanol

        # Test with show=False to avoid opening a window during tests
        fig = plot_assignment(smiles1, smiles2, assignment, show=False)

        # Check that the result is a matplotlib Figure
        assert isinstance(fig, plt.Figure)

    def test_plot_assignment_save(self, simple_molecules, assignment):
        """Test that plot_assignment can save the figure to a file."""
        mol1, mol2 = simple_molecules

        # Create a temporary file to save the figure to
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            # Test with show=False to avoid opening a window during tests
            plot_assignment(mol1, mol2, assignment, save_path=tmp_path, show=False)

            # Check that the file exists and has a non-zero size
            assert os.path.exists(tmp_path)
            assert os.path.getsize(tmp_path) > 0
        finally:
            # Clean up the temporary file
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
