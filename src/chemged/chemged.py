"""handles calculation of approximate chemical GEDs"""

import itertools
from typing import Literal, Sequence, overload

import numpy as np
import numpy.typing as npt
from networkx import Graph
from networkx.classes.coreviews import AtlasView
from scipy.optimize import linear_sum_assignment
from tqdm import tqdm

from .chem_utils import Molable, mol_to_nx, to_mol
from .cost import ChemicalGEDCostMatrix


MAX_INT32 = 2147483647


class ApproximateChemicalGED:
    r"""
    Calculator for approximate graph edit distance (GED) between two chemical graphs.

    Parameters
    ----------
    cost_matrix : ChemicalGEDCostMatrix
        the cost matrix that defines the costs of node and edge edit operations
        on the chemical graphs.

    Notes
    -----
    This implements a bipartite graph matching algorithm to approximate GED[1]_.

    The algorithm requires a way to compute the cost matrix of node and edge edit operations
    between two graphs. The cost matrix between two graphs is defined square matrix of size
    number of nodes in graph 1 + number of nodes in graph 2.

    .. math::

    \\mathbf{C} =
    \\begin{bmatrix}
    c_{1,1} & c_{1,2} & \\cdots & c_{1,m} & c_{1,\\epsilon} & \\infty & \\cdots & \\infty \\\\
    c_{2,1} & c_{2,2} & \\cdots & c_{2,m} & \\infty & c_{2,\\epsilon} & \\ddots & \\vdots \\\\
    \\vdots & \\vdots & \\ddots & \\vdots & \\vdots & \\ddots & \\ddots & \\infty \\\\
    c_{n,1} & c_{n,2} & \\cdots & c_{n,m} & \\infty & \\cdots & \\infty & c_{n,\\epsilon}\\\\
    c_{\\epsilon, 1} & \\infty & \\cdots & \\infty & 0 & 0 & \\cdots & 0 \\\\
    \\infty & c_{\\epsilon, 2} & \\ddots & \\vdots & 0 & 0 & \\ddots & \\vdots \\\\
    \\vdots & \\ddots & \\ddots & \\infty & \\vdots & \\ddots & \\ddots & 0 \\\\
    \\infty & \\cdots & \\infty & c_{\\epsilon, m} & 0 & \\cdots & 0 & 0 \\\\
    \\end{bmatrix}

    where :math:`c_{n,m}` is the cost of substituting node :math:`n` in graph :math:`n1`
    with node :math:`m` in graph :math:`n2`, :math:`c_{n,\\epsilon}` is the cost of deleting
    a node from :math:`n1`, and :math:`c_{\\epsilon,m}` is the cost of inserting node
    :math:`m` from :math:`n2` into :math:`n1`.

    Nodes can only be deleted or inserted once, so only the diagonals are populated in the
    insertion (bottom left) and deletion (upper right) sections of the cost matrix.
    The cost of doing another insertion or deletion is set to infinity, as it is not allowed.
    Particularly, we set it to the largest integer value to avoid numerical issues with the
    implementation of Munkres algorithm used by SciPy.
    The bottom right is always 0, as it represents the cost of substituting a missing node with
    another missing node, which has no cost.

    Then the Munkres algorithm is applied to find the optimal assignment of nodes and edges
    between the two graphs that minimizes the total cost of the edit operations.
    These costs are then summed to approximate the graph edit distance

    References
    ----------
    .. [1] Neuhaus, C., & Bunke, H. (2007). A graph edit distance
           based on the maximal common subgraph with application to
           structural document image recognition. In Document Analysis
           Systems VII (pp. 188-199). IEEE.
    """

    def __init__(self, cost_matrix: ChemicalGEDCostMatrix):
        self.cost_matrix = cost_matrix

    def _edge_cost_matrix(self, n1: AtlasView, n2: AtlasView):
        """
        Given two nodes, generate the cost matrix for edge edit operations.

        Parameters
        ----------
        n1: AtlasView
            the adj info the node from the first graph
        n2: AtlasView
            the adj info the node from the second graph

        Returns
        -------
        np.ndarray
            the cost matrix for edge edit operations between the two nodes
            shape will be (degree(n1) + degree(n2), degree(n1) + degree(n2))
        """
        cost_matrix = np.zeros([len(n1) + len(n2), len(n1) + len(n2)])

        # Insertion
        cost_matrix[len(n1) :, 0 : len(n2)] = MAX_INT32
        np.fill_diagonal(
            cost_matrix[len(n1) :, 0 : len(n2)], self.cost_matrix.get_edge_insertion_costs(n1)
        )

        # Deletion
        cost_matrix[0 : len(n1), len(n2) :] = MAX_INT32
        np.fill_diagonal(
            cost_matrix[0 : len(n1), len(n2) :], self.cost_matrix.get_edge_deletion_costs(n2)
        )

        # Substitution
        cost_matrix[0 : len(n1), 0 : len(n2)] = self.cost_matrix.get_edge_substitution_costs(
            n1, n2
        )

        return cost_matrix

    def _edge_edit_cost(self, n1: AtlasView, n2: AtlasView):
        """
        Compute the approximate graph edit cost for edges between two nodes.

        This is the part of the algorithm that attempts to match the local structures
        to avoid a much more expensive global edit distance calculation.

        Parameters
        ----------
        n1: AtlasView
            the adj info the node from the first graph
        n2: AtlasView
            the adj info the node from the second graph

        Returns
        -------
        float
            the approximate graph edit cost for edges between the two nodes
        """
        # Compute cost matrix
        cost_matrix = self._edge_cost_matrix(n1, n2)

        # Munkres algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Edge edit distance between the two nodes
        dist = cost_matrix[row_ind, col_ind].sum()

        return dist

    def _get_full_cost_matrix(self, g1: Graph, g2: Graph):
        """
        Compute the full cost matrix for the two graphs, including node and edge edit operations.

        This matrix is what is used to compute the approximate graph edit distance with
        the Munkres algorithm.

        Parameters
        ----------
        g1: Graph
            the first graph
        g2: Graph
            the second graph

        Returns
        -------
        np.ndarray
            the full cost matrix for the two graphs, shape will be
            (len(g1) + len(g2), len(g1) + len(g2)) where len(g) is
            the number of nodes in the graph.
        """
        cost_matrix = np.zeros([len(g1) + len(g2), len(g1) + len(g2)])

        # Insertion
        cost_matrix[len(g1) :, 0 : len(g2)] = MAX_INT32
        np.fill_diagonal(
            cost_matrix[len(g1) :, 0 : len(g2)],
            self.cost_matrix.get_node_insertion_costs(g1)
            + [self.cost_matrix.get_edge_insertion_costs(g1[i]).sum() for i in range(len(g1))],
        )

        # Deletion
        cost_matrix[0 : len(g1), len(g2) :] = MAX_INT32
        np.fill_diagonal(
            cost_matrix[0 : len(g1), len(g2) :],
            self.cost_matrix.get_node_deletion_costs(g2)
            + [self.cost_matrix.get_edge_deletion_costs(g2[i]).sum() for i in range(len(g2))],
        )

        # Substitution
        node_dist = self.cost_matrix.get_node_substitution_costs(g1, g2)
        # update the substitution costs with the edge edit scores
        for i, n in enumerate(g1.nodes()):
            for j, m in enumerate(g2.nodes()):
                node_dist[i, j] += self._edge_edit_cost(g1[n], g2[m])
        cost_matrix[0 : len(g1), 0 : len(g2)] = node_dist

        return cost_matrix

    @overload
    def _ged(
        self, g1: Graph, g2: Graph, return_assignment: Literal[False]
    ) -> tuple[float, None]: ...

    @overload
    def _ged(
        self, g1: Graph, g2: Graph, return_assignment: Literal[True]
    ) -> tuple[float, tuple[npt.NDArray, npt.NDArray]]: ...

    @overload
    def _ged(
        self, g1: Graph, g2: Graph, return_assignment: bool
    ) -> tuple[float, tuple[npt.NDArray, npt.NDArray] | None]: ...

    def _ged(
        self, g1: Graph, g2: Graph, return_assignment: bool = False
    ) -> tuple[float, tuple[npt.NDArray, npt.NDArray] | None]:
        """
        Compute the approximate graph edit distance between two chemicals

        Parameters
        ----------
        g1: Graph
            the first chemical graph
        g2: Graph
            the second chemical graph

        Returns
        -------
        distance, assignment: float, tuple[npt.ndarray, npt.ndarray] | None
            'distance' is the approximate graph edit distance between the two graphs,
            'assignment' is a tuple of two numpy arrays defining the assignment
            between the nodes of the two graphs
        """
        # Compute cost matrix
        cost_matrix = self._get_full_cost_matrix(g1, g2)

        # Munkres algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Graph edit distance
        dist = cost_matrix[row_ind, col_ind].sum()

        if return_assignment:
            not_assign = np.invert((row_ind >= len(g1)) * (col_ind >= len(g2)))
            return dist, (row_ind[not_assign], col_ind[not_assign])
        else:
            return dist, None

    @overload
    def compute_ged(
        self, chemical1: Molable, chemical2: Molable, return_assignment: Literal[False]
    ) -> float: ...

    @overload
    def compute_ged(
        self, chemical1: Molable, chemical2: Molable, return_assignment: Literal[True]
    ) -> tuple[float, tuple[npt.NDArray, npt.NDArray]]: ...

    @overload
    def compute_ged(
        self, chemical1: Molable, chemical2: Molable, return_assignment: bool
    ) -> float | tuple[float, tuple[npt.NDArray, npt.NDArray]]: ...

    def compute_ged(
        self, chemical1: Molable, chemical2: Molable, return_assignment: bool = False
    ) -> float | tuple[float, tuple[npt.NDArray, npt.NDArray]]:
        """
        Compute the approximate graph edit distance between two chemicals

        Will convert the chemicals to rdkit.Mol objects if they are not already,
        and then covert them into NetworkX Graph objects. See `mol_to_nx` for more details.

        NOTE: it is inefficient to call this function multiple times with the same chemicals,
        as the chemicals will be converted to NetworkX Graph objects each time. If you have
        a set of chemcials that you want to get a pairwise distance for, use the `pdist``
        or `cdist` functions instead.

        Parameters
        ----------
        chemical1: Molable (str or rdkit.Mol)
            the first chemical to compare, can be a SMILES string or rdkit.Mol object
        chemical2: Molable (str or rdkit.Mol)
            the second chemical to compare, can be a SMILES string or rdkit.Mol object
        return_assignment: bool, default=False
            whether to return the assignment of nodes between the two graphs.
            If True, will return a tuple of two numpy arrays defining the assignment
            between the nodes of the two graphs. If False, will only return the distance.

        Returns
        -------
        distance, assignment: float, tuple[float, tuple[npt.ndarray, npt.ndarray]]
            'distance' is the approximate graph edit distance between the two chemicals,
            'assignment' is a tuple of two numpy arrays defining the assignment
            between the nodes of the two graphs if `return_assignment` is True.
        """
        mol1 = to_mol(chemical1, fail_on_error=True)
        mol2 = to_mol(chemical2, fail_on_error=True)

        g1: Graph = mol_to_nx(mol1)
        g2: Graph = mol_to_nx(mol2)

        dist, assignment = self._ged(g1, g2, return_assignment=return_assignment)
        if return_assignment:
            if assignment is None:
                raise ValueError("Assignment is None, this should not happen.")
            else:
                return dist, assignment
        else:
            return dist

    def pdist(self, chemicals: Sequence[Molable], use_tqdm: bool = False) -> npt.NDArray:
        """
        Given a list of chemicals, get the approximate graph edit distance between chemicals

        Notes
        -----
        This will return distances in vector-form distance vector.
        Yu can use `scipy.spatial.distance.squareform` to convert this to a square distance matrix.

        Parameters
        ----------
        chemicals: Sequence[Molable]
            a sequence of chemicals to compute the pairwise distances for.
            Each chemical can be a SMILES string or rdkit.Mol object.
        use_tqdm: bool, default=False
            enable a tqdm progress bar for the pairwise distance calculation.

        Returns
        -------
        np.ndarray
            a condensed distance vector containing the pairwise distances between the chemicals.
            The length of the vector will be n * (n - 1) / 2, where n is the number of chemicals.
        """
        mols = [
            to_mol(c, fail_on_error=True)
            for c in tqdm(
                chemicals, disable=not use_tqdm, desc="Converting chemicals to Mol objects"
            )
        ]
        graphs = [
            mol_to_nx(mol)
            for mol in tqdm(
                mols, disable=not use_tqdm, desc="Converting Mol objects to NetworkX graphs"
            )
        ]

        dists = []
        total_calcs = len(graphs) * (len(graphs) - 1) // 2
        for g1, g2 in tqdm(
            itertools.combinations(graphs, 2),
            total=total_calcs,
            disable=not use_tqdm,
            desc="Computing pairwise Approx. GED distances between chemicals",
        ):
            dists.append(self._ged(g1, g2, return_assignment=False)[0])

        return np.array(dists, dtype=np.float64)

    def cdist(
        self, chemicals1: Sequence[Molable], chemicals2: Sequence[Molable], use_tqdm: bool = False
    ) -> npt.NDArray:
        """
        Given two sequences of chemicals, compute the approximate graph edit distance between them

        Parameters
        ----------
        chemicals1: Sequence[Molable]
            a sequence of chemicals to compute the pairwise distances for.
            Each chemical can be a SMILES string or rdkit.Mol object.
        chemicals2: Sequence[Molable]
            a second sequence of chemicals to compute the pairwise distances for.
            Each chemical can be a SMILES string or rdkit.Mol object.
        use_tqdm: bool, default=False
            enable a tqdm progress bar for the pairwise distance calculation.

        Returns
        -------
        np.ndarray
            a condensed distance vector containing the pairwise distances between the chemicals.
            The length of the vector will be n * m, where n is the number of chemicals in
            `chemicals1` and m is the number of chemicals in `chemicals2`.
        """
        mols1 = [
            to_mol(c, fail_on_error=True)
            for c in tqdm(
                chemicals1,
                disable=not use_tqdm,
                desc="Converting first set of chemicals to Mol objects",
            )
        ]
        mols2 = [
            to_mol(c, fail_on_error=True)
            for c in tqdm(
                chemicals2,
                disable=not use_tqdm,
                desc="Converting second set of chemicals to Mol objects",
            )
        ]

        graphs1 = [
            mol_to_nx(mol)
            for mol in tqdm(
                mols1,
                disable=not use_tqdm,
                desc="Converting first set of Mol objects to NetworkX graphs",
            )
        ]
        graphs2 = [
            mol_to_nx(mol)
            for mol in tqdm(
                mols2,
                disable=not use_tqdm,
                desc="Converting second set of Mol objects to NetworkX graphs",
            )
        ]

        dists = []
        total_calcs = len(graphs1) * len(graphs2)
        for g1, g2 in tqdm(
            itertools.product(graphs1, graphs2),
            total=total_calcs,
            disable=not use_tqdm,
            desc="Computing pairwise Approx. GED distances between chemicals",
        ):
            dists.append(self._ged(g1, g2, return_assignment=False)[0])

        return np.array(dists, dtype=np.float64)
