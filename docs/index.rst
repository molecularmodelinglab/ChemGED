=======
ChemGED
=======

ChemGED is a Python package for enabling the appoximate graph edit distance (GED) computation between
chemicals. Normally, GED is a NP-hard problem, but ChemGED uses heuristics to approximate the GED in
a reasonable time.

.. code-block:: python

    from chemged import ApproximateChemicalGED, UniformElementCostMatrix

    ged_calc = ApproximateChemicalGED(cost_matrix=UniformElementCostMatrix())

    # Compute the approximate GED
    ged = ged_calc.compute_ged("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "CC(=O)Nc1ccc(O)cc1")
    print(f"Approximate GED: {ged}")
    >>> Approximate GED: 12.0


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   chemical_graph_attributes
   custom_cost_matrix
   api
   contributing

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
