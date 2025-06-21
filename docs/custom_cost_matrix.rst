.. _custom_cost_matrix:

====================
Custom Cost Matrices
====================

The cost matrix is a crucial component of the Approximate Graph Edit Distance (GED) calculation used in ChemGED.

This matrix takes form:

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
:math:`m` from :math:`n2` into :math:`n1`. Note that in this cost matrix form, the
cost of substituting a node takes into account the cost of the edge operations as well.
The edge cost is also defined using the same matrix structure, except for the edge
between two nodes. A more detailed explanation of the cost matrix can be found in the
paper by [1]_.

This is alot to handle, so ChemGED provides a abstract class ``ChemicalGEDCostMatrix``
that defines the interface for cost matrices definition. Rather than defining how to make
the full cost matrix in one go. Instead you need to define the follow cost functions:

- Node substitution: Generates a cost matrix for replacing a node in one graph with a node in another graph
  (shape: Num Nodes in Graph 1 x Num Nodes in Graph 2)
- Node insertion: Generates a cost vector for inserting a node from from the graph (shape: Num Nodes in Graph)
- Node deletion: Generates a cost vector for deleting a node from from the graph (shape: Num Nodes in Graph)
- Edge substitution: Generates a cost matrix for replacing all the edges in one node with edges in another node
  (shape: Num Edges in node 1 x Num Edges in node 2)
- Edge insertion: Generates a cost vector for inserting a edge from from the node (shape: Num edges in Graph)
- Edge deletion: Generates a cost vector for deleting a edge from from the graph (shape: Num edge in Graph)

With these functions defined, ChemGED will take care of the rest to construct the full cost matrix and run the
approximate GED algorithm.

By default, ChemGED will use the ``UniformElementCostMatrix`` which defines a uniform cost for all operations.
This is a good starting point, but you allows limited control over atomic properties and bond types. Swapping
a Cl with a Br will have a bigger impact of the chemical than swapping a C with a N, for example. If this is
something you want to take into account, you can create your own cost matrix by subclassing
``ChemicalGEDCostMatrix`` and implementing the required methods taking into account the chemical properties of
the nodes and edges.

Creating a Custom Cost Matrix
=============================

To create a custom cost matrix, you need to subclass ``ChemicalGEDCostMatrix`` and implement the required methods:

.. code-block:: python

    from chemged import ChemicalGEDCostMatrix
    import networkx as nx
    import numpy as np

    class MyCustomCostMatrix(ChemicalGEDCostMatrix):
        def __init__(self, custom_param1, custom_param2):
            # Initialize your custom parameters if needed
            self.custom_param1 = custom_param1
            self.custom_param2 = custom_param2

        def get_node_substitution_costs(self, g1: nx.Graph, g2: nx.Graph) -> np.ndarray:
            # Implement your custom node substitution cost logic
            # Return a matrix of shape (len(g1.nodes), len(g2.nodes))
            # where each element (i, j) is the cost of substituting node i in g1 with node j in g2
            pass

        def get_node_insertion_costs(self, g: nx.Graph) -> np.ndarray:
            # Implement your custom node insertion cost logic
            # Return a vector of shape (len(g.nodes),)
            # where each element i is the cost of inserting node i
            pass

        def get_node_deletion_costs(self, g: nx.Graph) -> np.ndarray:
            # Implement your custom node deletion cost logic
            # Return a vector of shape (len(g.nodes),)
            # where each element i is the cost of deleting node i
            pass

        def get_edge_substitution_costs(self, n1: nx.classes.coreviews.AtlasView, n2: nx.classes.coreviews.AtlasView) -> np.ndarray:
            # Implement your custom edge substitution cost logic
            # Return a matrix of shape (len(n1), len(n2))
            # where each element (i, j) is the cost of substituting edge i in n1 with edge j in n2
            pass

        def get_edge_insertion_costs(self, n: nx.classes.coreviews.AtlasView) -> np.ndarray:
            # Implement your custom edge insertion cost logic
            # Return a vector of shape (len(n),)
            # where each element i is the cost of inserting edge i
            pass

        def get_edge_deletion_costs(self, n: nx.classes.coreviews.AtlasView) -> np.ndarray:
            # Implement your custom edge deletion cost logic
            # Return a vector of shape (len(n),)
            # where each element i is the cost of deleting edge i
            pass

Example: Custom Cost Matrix to Ignore Halogen Swaps
===================================================

Here's an example of how to define the a custom cost matrix that modifies the
UniformElementCostMatrix to ignore swaps between halogens (F, Cl, Br, I):

.. code-block:: python

    from chemged import UniformElementCostMatrix
    import networkx as nx

    class HalogenSwapCostMatrix(UniformElementCostMatrix):
        def get_node_substitution_costs(self, g1: nx.Graph, g2: nx.Graph) -> np.ndarray:
            v1 = list(nx.get_node_attributes(g1, "atomic_num").values())
            v2 = list(nx.get_node_attributes(g2, "atomic_num").values())

            v1[np.where(np.isin(v1, [9, 17, 35, 53]))] = -1  # alias for F, Cl, Br, I
            v2[np.where(np.isin(v2, [9, 17, 35, 53]))] = -1  # alias for F, Cl, Br, I

            # generate the cost matrix
            return ~np.equal.outer(v1, v2) * self.node_sub_cost

Now the chemicals "CCCCCCl" and "CCCCCBr" will have a GED of 0.
This can be a very helpful tool if trying to build more chemically aware
similarity calculations.

.. note::
    atoms in the networkx graph have a handful of attributes other than "atomic_num", see the
    :ref:`chemical graph attributes <chem-graph-attributes>` for more information
    on what you can use (or how to add new ones).

Using Your Custom Cost Matrix
=============================

Once you've implemented your custom cost matrix, you can use it with the ApproximateChemicalGED class:

.. code-block:: python

    from chemged import ApproximateChemicalGED
    from my_module import HalogenSwapCostMatrix

    # Create an instance of your custom cost matrix
    custom_cost_matrix = HalogenSwapCostMatrix()

    # Create an instance of ApproximateChemicalGED with your custom cost matrix
    ged_calc = ApproximateChemicalGED(cost_matrix=custom_cost_matrix)

    # Compute the approximate GED
    ged = ged_calc.compute_ged("CCCCCCl", "CCCCCBr")
    print(f"Approximate GED: {ged}")
    >>> Approximate GED: 0.0

Tips for Implementing Custom Cost Matrices
==========================================
While the approximate GED algorithm is far faster than a brute for GED calculation, it still can be computationally expensive,
especially if generating the cost matrix is not efficient. Try to limit the number of operations in the cost matrix generation
functions to only what is necessary for you case, and keep it as effiecent as you can.

References
==========
.. [1] Riesen, Kaspar, and Horst Bunke. "Approximate graph edit distance computation by means of bipartite graph matching." Image and Vision computing 27.7 (2009): 950-959
