# ChemGED

ChemGED is a Python package for enabling the appoximate graph edit distance (GED) computation between 
chemicals. Normally, GED is a NP-hard problem, but ChemGED uses heuristics to approximate the GED in 
a reasonable time.

## Installation
You can install ChemGED using pip:

```bash
pip install chemged
```

## Usage
To use ChemGED, you just create an ``ApproximateChemicalGED`` object, and then call the 
``compute_ged`` method with two chemical structures. They can be SMILES *or* 
[RDKit](https://www.rdkit.org/docs/index.html) Mol objects.

```python
from chemged import ApproximateChemicalGED
ged_calc = ApproximateChemicalGED()

chemical1 = "CCO"  # SMILES of the first chemical
chemical2 = "CCN"  # SMILES of the second chemical

# you can use SMILES strings directly
approx_ged = ged_calc.compute_ged(chemical1, chemical2)

# you can also use RDKit Mol objects
from rdkit.Chem import MolFromSmiles
approx_ged = ged_calc.compute_ged(MolFromSmiles(chemical1), MolFromSmiles(chemical2))
```
ChemGED also implements ``pdist`` and ``cdist`` functions to compute pairwise distances between
sets of chemicals. These will return as Numpy arrays.

> [!NOTE]  
> ``pdist`` will return the vector-form distance vector, while ``cdist`` will return a
> square-form distance matrix. ``scipy.spatial.distance.squareform`` can be used to convert
> the vector-form distance vector to a square-form distance matrix.

```python

```


## Implementation

The approach used here uses bipartite graph matching[[1]](#1), and most of its implementation in python 
is based off scripts from https://github.com/priba/aproximated_ged/tree/master.
ChemGED uses [RDKit](https://www.rdkit.org/docs/index.html) to handle chemicals inside python.


## References
<a id="1">[1]</a> 
Riesen, Kaspar, and Horst Bunke. 
"Approximate graph edit distance computation by means of bipartite graph matching." 
Image and Vision computing 27.7 (2009): 950-959