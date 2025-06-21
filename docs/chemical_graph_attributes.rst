.. _chem-graph-attributes:

=========================
Chemical Graph Attributes
=========================

ChemGED uses the NetworkX library to represent chemical graphs. Each node in the graph represents an atom,
and each edge represents a bond between atoms. The attributes of these nodes and edges are crucial for
calculating the graph edit distance (GED) accurately.

By default, the ``mol_to_nx`` function in ChemGED will convert a chemical into a NetworkX graph with the following
atom attributes stored in the node:

- **atomic_number**: The atomic number of the atom (e.g., 6 for Carbon, 8 for Oxygen).
- **degree**: The degree of the atom in the graph (number of bonds).
- **formal_charge**: The formal charge of the atom.
- **hybridization**: The hybridization state of the atom (e.g., 'sp3', 'sp2').
- **is_aromatic**: A boolean indicating if the atom is part of an aromatic ring.
- **is_in_ring**: A boolean indicating if the atom is part of a ring structure.
- **hydrogen**: The number of hydrogen atoms attached to the atom.

Chemical bond attributes are stored in the graph. The edge attributes include:
- **bond_type**: The type of bond (e.g., 'single', 'double', 'triple', 'aromatic').

While these are the current attributes stored in the graph, you can extend or modify these attributes as
needed for your specific use case.
