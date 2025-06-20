"""For handling the cost matrix for chemical graph edit distance."""

import abc

import numpy as np
import numpy.typing as npt

import networkx as nx
from networkx import Graph
from networkx.classes.coreviews import AtlasView


class ChemicalGEDCostMatrix(abc.ABC):
    """
    Abstract base class for chemical graph edit distance cost matrix

    These Matrices are used to determine the cost of any edit operation.

    The costs of node substitution, insertion and deletion costs
    need to be defined. These should return a cost matrix (for substitution)
    or a cost vector (for insertion and deletion) of the nodes.
    This is done by the functions
    - `get_node_substitution_costs`
    - `get_node_insertion_costs`
    - `get_node_deletion_costs`

    Likewise, the cost functions for edge substitution, insertion and deletion
    are needed. These are a bit different. Rather than operating on the whole graph,
    they operate locally on the nodes. They should return a cost matrix (substitution)
    or a cost vector (insertion and deletion) of the edges *for the two given nodes*
    only. So if we have two nodes n1 (from graph1) and n2 (from graph2), and n1 has
    a degree of 3 (3 bonds) and n2 has a degree of 2 (2 bonds), then the substitution
    cost matrix should be of size 3x2, where the element at position (i, j) is the
    substitution cost of the i-th bond of n1 with the j-th bond of n2. Likewise,
    the deletion/insertion cost vector should be the same size as the respective node's degree.

    - `get_edge_substitution_costs`
    - `get_edge_insertion_costs`
    - `get_edge_deletion_costs`
    """

    @abc.abstractmethod
    def get_node_substitution_costs(self, g1: Graph, g2: Graph):
        """
        Get the node substitution costs as a matrix between two ordered graphs.

        This method should return a cost matrix where the element at position (i, j)
        represents the cost of substituting node i in n1 with node j in n2.

        The node order in these graphs *MUST* stay static over the object lifetime
        (at least during GED calculation),
        otherwise the generated cost matrix may become invalid during calculation.

        Parameters
        ----------
        g1 : Graph
            The first chemical graph.
        g2 : Graph
            The second chemical graph.

        Returns
        -------
        np.ndarray
            A 2D cost matrix for node substitutions.
            shape should be (len(n1.nodes), len(n2.nodes))
        """
        raise NotImplementedError

    @abc.abstractmethod
    def get_edge_substitution_costs(self, n1: AtlasView, n2: AtlasView):
        """
        Get the edge substitution costs as a matrix between two nodes.

        This method should return a cost matrix where the element at position (i, j)
        represents the cost of substituting edge i in n1 with edge j in n2.

        The edge order in these graphs *MUST* stay static over the object lifetime
        (at least during GED calculation),
        otherwise the generated cost matrix may become invalid during calculation.

        Notes
        -----
        You can get the adjacency list of a node in a graph from network x using
        graph.adj[n] OR using g[n] where n is the node key/index. These are called
        `AtlasViews` in NetworkX and can be treated as dictionaries for the most part.

        Parameters
        ----------
        n1 : AtlasView
            The adjacency list of the first chemical graph node.
        n2 : AtlasView
            The adjacency list of the first chemical graph node.

        Returns
        -------
        np.ndarray
            A 2D cost matrix for edge substitutions.
            shape should be (degree(n1), degree(n2))
        """
        raise NotImplementedError

    def get_node_insertion_costs(self, g: Graph):
        """
        Get the node insertion costs as a vector for a graph.

        This method should return a cost vector where the element at position i
        represents the cost of inserting node i in g.

        Parameters
        ----------
        g : Graph
            The chemical graph.

        Returns
        -------
        np.ndarray
            A 1D cost vector for node insertions.
            shape should be (len(g.nodes),)
        """
        raise NotImplementedError

    def get_node_deletion_costs(self, g: Graph):
        """
        Get the node deletion costs as a vector for a graph.

        This method should return a cost vector where the element at position i
        represents the cost of deleting node i in g.

        Parameters
        ----------
        g : Graph
            The chemical graph.

        Returns
        -------
        np.ndarray
            A 1D cost vector for node deletions.
            shape should be (len(g.nodes),)
        """
        raise NotImplementedError

    def get_edge_insertion_costs(self, n: AtlasView):
        """
        Get the edge insertion costs as a vector for a node in a graph.

        This method should return a cost vector where the element at position i
        represents the cost of inserting edge i in the node n.

        Parameters
        ----------
        n : AtlasView
            The adjacency list of the chemical graph node.

        Returns
        -------
        np.ndarray
            A 1D cost vector for edge insertions.
            shape should be (degree(n),)
        """
        raise NotImplementedError

    def get_edge_deletion_costs(self, n: AtlasView):
        """
        Get the edge deletion costs as a vector for a node in a graph.

        This method should return a cost vector where the element at position i
        represents the cost of deleting edge i in the node n.

        Parameters
        ----------
        n : AtlasView
            The adjacency list of the chemical graph node.

        Returns
        -------
        np.ndarray
            A 1D cost vector for edge deletions.
            shape should be (degree(n),)
        """
        raise NotImplementedError


class UniformElementCostMatrix(ChemicalGEDCostMatrix):
    """
    A uniform cost matrix for chemical graph edit distance.

    This cost matrix assumes all substitutions, insertions, and deletions of atoms
    are worth a cost of 1, not matter which type the atom is *or* what type the
    atom was swapped into. The same goes for bonds, where all bond substitutions,
    insertions, and deletions are worth a cost of 1, regardless of the bond type.

    This is good for if you want to measure GED in a more "pure" sense, where
    we are just asking how much the graph representing a chemical has changed.
    If you want to account for chemical properties, you should write a more
    advanced cost matrix that takes into account the atom and bond types.
    For example, minimize cost for O to N substitutions if they are both
    Hydrogen bond donors, but increase the cost for O to C substitutions.

    Parameters
    ----------
    node_sub_cost : float, default=1.0
        The uniform cost of substituting a node (atom).
    node_ins_cost : float, default=1.0
        The uniform cost of inserting a node (atom).
    node_del_cost : float, default=1.0
        The uniform cost of deleting a node (atom).
    edge_sub_cost : float, default=1.0
        The uniform cost of substituting an edge (bond).
    edge_ins_cost : float, default=1.0
        The uniform cost of inserting an edge (bond).
    edge_del_cost : float, default=1.0
        The uniform cost of deleting an edge (bond).
    """
    def __init__(
            self,
            node_sub_cost=1.0,
            node_ins_cost=1.0,
            node_del_cost=1.0,
            edge_sub_cost=1.0,
            edge_ins_cost=1.0,
            edge_del_cost=1.0
    ):
        """initialize an object"""
        self.node_sub_cost = node_sub_cost
        self.node_ins_cost = node_ins_cost
        self.node_del_cost = node_del_cost

        self.edge_sub_cost = edge_sub_cost
        self.edge_ins_cost = edge_ins_cost
        self.edge_del_cost = edge_del_cost

    def get_node_substitution_costs(self, g1: Graph, g2: Graph) -> npt.NDArray:
        """
        Get the node substitution costs as a matrix between two ordered graphs.

        This method returns a cost matrix where the element at position (i, j)
        represents the cost of substituting node i in n1 with node j in n2.

        It uses only the element type to determine if two nodes are the same,
        if they are the cost is 0, otherwise it is `self.node_sub_cost`.

        Parameters
        ----------
        g1 : Graph
            The first chemical graph.
        g2 : Graph
            The second chemical graph.

        Returns
        -------
        np.ndarray
            A 2D cost matrix for node substitutions.
            shape should be (len(n1.nodes), len(n2.nodes))
        """
        # all the atomic numbers of the nodes in the graphs
        v1 = list(nx.get_node_attributes(g1, "atomic_num").values())
        v2 = list(nx.get_node_attributes(g2, "atomic_num").values())

        # generate the cost matrix
        return np.equal.outer(v1, v2) * self.node_sub_cost

    def get_node_insertion_costs(self, g: Graph) -> npt.NDArray:
        """
        Get the node insertion costs as a vector for a graph.

        Will return a vector of the same length as the number of nodes in the graph,
        filled with the value in `self.node_ins_cost.

        Parameters
        ----------
        g: Graph
            the chemical graph to get the insertion costs for

        Returns
        -------
        np.ndarray
            A 1D cost vector for node insertions.
            shape should be (len(g.nodes),)
        """
        return np.full(len(g), self.node_ins_cost)

    def get_node_deletion_costs(self, g: Graph) -> npt.NDArray:
        """
        Get the node deletion costs as a vector for a graph.

        Will return a vector of the same length as the number of nodes in the graph,
        filled with the value in `self.node_del_cost.

        Parameters
        ----------
        g: Graph
            the chemical graph to get the deletion costs for

        Returns
        -------
        np.ndarray
            A 1D cost vector for node deletions.
            shape should be (len(g.nodes),)
        """
        return np.full(len(g), self.node_del_cost)

    def get_edge_substitution_costs(self, n1: AtlasView, n2: AtlasView) -> npt.NDArray:
        """
        Get the edge substitution costs as a matrix between two nodes.

        This method returns a cost matrix where the element at position (i, j)
        represents the cost of substituting edge i in n1 with edge j in n2.

        It uses only the bond type to determine if two edges are the same,
        if they are the cost is 0, otherwise it is `self.edge_sub_cost`.

        Parameters
        ----------
        n1 : AtlasView
            The adjacency list of the first chemical graph node.
        n2 : AtlasView
            The adjacency list of the second chemical graph node.

        Returns
        -------
        np.ndarray
            A 2D cost matrix for edge substitutions.
            shape should be (degree(n1), degree(n2))
        """
        # get the bond types of the edges in the nodes
        e1 = [e["bond_type"] for e in n1.values()]
        e2 = [e["bond_type"] for e in n2.values()]

        # generate the cost matrix
        return np.equal.outer(e1, e2) * self.edge_sub_cost

    def get_edge_insertion_costs(self, n: AtlasView) -> npt.NDArray:
        """
        Get the edge insertion costs as a vector for a node in a graph.

        Will return a vector of the same length as the number of edges in the node,
        filled with the value in `self.edge_ins_cost`.

        Parameters
        ----------
        n : AtlasView
            The adjacency list of the chemical graph node.

        Returns
        -------
        np.ndarray
            A 1D cost vector for edge insertions.
            shape should be (degree(n),)
        """
        return np.full(len(n), self.edge_ins_cost)

    def get_edge_deletion_costs(self, n: AtlasView) -> npt.NDArray:
        """
        Get the edge deletion costs as a vector for a node in a graph.

        Will return a vector of the same length as the number of edges in the node,
        filled with the value in `self.edge_del_cost`.

        Parameters
        ----------
        n : AtlasView
            The adjacency list of the chemical graph node.

        Returns
        -------
        np.ndarray
            A 1D cost vector for edge deletions.
            shape should be (degree(n),)
        """
        return np.full(len(n), self.edge_del_cost)
