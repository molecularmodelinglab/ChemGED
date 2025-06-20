try:
    import matplotlib.pyplot as plt
except ImportError as e:
    raise ImportError("matplotlib is required for plotting. Please install it using 'pip install matplotlib'.") from e

from os import PathLike
from typing import Optional

import networkx as nx
import numpy as np

from .chem_utils import Molable, to_mol, mol_to_nx


def plot_assignment(chemical1:Molable, chemical2: Molable, assignment, save_path: Optional[PathLike]=None, show: bool=True) -> "plt.Figure":
    """
    Plot the assignment of two chemical structures, mapping nodes from one to the other.

    This might be useful to visualize if the approximate GED calculator is making good
    node mapping. It could be hard to interpret the results if there are lots of atoms
    though.

    Parameters
    ----------
    chemical1: Molable
        the first chemical structure, can be a SMILES string or an RDKit Mol object
        MUST be the same chemical as the one used to generate the assignment
    chemical2: Molable
        the second chemical structure, can be a SMILES string or an RDKit Mol object
        MUST be the same chemical as the one used to generate the assignment
    assignment: tuple[np.ndarray, np.ndarray]
        the assignment of two chemical structures, mapping nodes from one to the other
        returned by the approximate GED calculator if "return_assignment" is True
    save_path: Optional[PathLike], default=None
        path to save the plot to, if None, the plot will not be saved
    show: bool, default=True
        show the plot after generating it, if False, the plot will not be shown

    Returns
    -------
    Figure
        the matplotlib Figure object containing the plot of the assignment
    """
    g1 =  mol_to_nx(to_mol(chemical1, fail_on_error=True))
    g2 =  mol_to_nx(to_mol(chemical2, fail_on_error=True))

    g_comp = nx.disjoint_union(g1, g2)

    position = {k: v['coord'] for k, v in g_comp.node.items()}

    center1 = np.mean(list(position.values())[0:len(g1)], axis=0)
    max_pos1 = np.max(np.abs(list(position.values())[0:len(g1)] - center1))

    center2 = np.mean(list(position.values())[len(g1):], axis=0)
    max_pos2 = np.max(np.abs(list(position.values())[len(g1):] - center2))

    keys_modify = list(position.keys())[len(g1):]
    for x in keys_modify:
        position[x][0] += 2 * (np.max([max_pos1, max_pos2]) + 0.5)

    ag1, ag2 = assignment

    edgelist = []
    nodelist_ins = []
    nodelist_del = []
    for i1, i2 in zip(ag1, ag2):
        if i1 < len(g1):
            if i2 < len(g2):
                # Substitution
                g_comp.add_edge(i1, i2 + len(g1))
                edgelist += [(i1, i2 + len(g1))]
            else:
                # Deletion
                nodelist_del += [i1]
        else:
            # Insertion
            nodelist_ins += [i2 + len(g1)]

    center = np.mean(list(position.values()), axis=0)
    max_pos = np.max(np.abs(np.array(position.values()) - center))

    fig = plt.figure()

    nx.draw_networkx_nodes(g_comp, position,
                           nodelist=[item for item in g_comp.nodes() if item not in nodelist_ins + nodelist_del],
                           node_color='black',
                           node_size=200)

    nx.draw_networkx_edges(g_comp, position,
                           edgelist=[item for item in g_comp.edges() if item not in edgelist]
                           )

    nx.draw_networkx_edges(g_comp, position,
                           edgelist=edgelist,
                           width=3, alpha=0.5, edge_color='b',
                           style='dashed')

    nx.draw_networkx_nodes(g_comp, position,
                           nodelist=nodelist_ins,
                           node_color='g',
                           node_size=500,
                           alpha=0.8)

    nx.draw_networkx_nodes(g_comp, position,
                           nodelist=nodelist_del,
                           node_color='r',
                           node_size=500,
                           alpha=0.8)

    plt.ylim([center[1] - max_pos - 0.5, center[1] + max_pos + 0.5])
    plt.xlim([center[0] - max_pos - 0.5, center[0] + max_pos + 0.5])
    plt.axis('off')

    if save_path is not None:
        plt.savefig(save_path)

    if show:
        plt.show()

    return fig
