"""
Manages hypergraphs
"""

from hypernetx import Hypergraph

from matplotlib.collections import PolyCollection

from hypernetx import Hypergraph
from hypernetx.drawing.util import (
    get_frozenset_label,
    get_collapsed_size,
    get_set_layering,
    inflate_kwargs,
)
from hypernetx.drawing.rubber_band import (
    layout_node_link,
    get_default_radius,
    layout_hyper_edges,
    draw_hyper_edge_labels,
    draw_hyper_labels,
)
import hypernetx as hnx
import hypernetx.reports.descriptive_stats as hnx_stats
import matplotlib.colors as mc
from matplotlib import colormaps as cm
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.spatial import ConvexHull

from microbial_hypergraphs.population import Population
from microbial_hypergraphs.correlation import Correlation
from microbial_hypergraphs.hypercorrelation import HyperCorrelation
from microbial_hypergraphs.core import LOGGER
from microbial_hypergraphs.core.data import get_taxonomy, get_otu_samples

###########################
### HYPERNETX OVERRIDES ###
###########################
# The hypernetx package lacks some of the coloring functionality that we would like to do.
# So I've brought in some portions of their package here and overwritten its functionality
# to include what we want to do. The method docstrings are from hypernetx.

N_CONTROL_POINTS = 24
theta = np.linspace(0, 2 * np.pi, N_CONTROL_POINTS + 1)[:-1]
cp = np.vstack((np.cos(theta), np.sin(theta))).T


def _draw_hyper_nodes(
    H, pos, node_color_map=cm["Greens"], node_radius={}, r0=None, ax=None, **kwargs
):
    """
    Draws a circle for each node in H.

    --- CUSTOM DOCS ---
    The

    --- HYPERNETX PACKAGE DOCS ---

    The position of each node is specified by the a dictionary/list-like, pos,
    where pos[v] is the xy-coordinate for the vertex. The radius of each node
    can be specified as a dictionary where node_radius[v] is the radius. If a
    node is missing from this dictionary, or the node_radius is not specified at
    all, a sensible default radius is chosen based on distances between nodes
    given by pos.

    Parameters
    ----------
    H: Hypergraph
        the entity to be drawn
    pos: dict
        mapping of node and edge positions to R^2
    node_radius: dict
        mapping of node to R^1 (radius of each node)
    r0: float
        minimum distance that concentric rings start from the node position
    ax: Axis
        matplotlib axis on which the plot is rendered
    kwargs: dict
        keyword arguments, e.g., linewidth, facecolors, are passed through to the PolyCollection constructor

    Returns
    -------
    PolyCollection
        a Matplotlib PolyCollection that can be further styled
    """

    taxonomy = get_taxonomy()
    otu_overall_presence = get_otu_samples().T.sum(axis=1)

    ax = ax or plt.gca()

    r0 = r0 or get_default_radius(H, pos)

    points = []
    colors = []
    for v in H.nodes:
        colors.append(otu_overall_presence.loc[v])
        match taxonomy.loc[v]["Domain"]:
            case "Archaea":
                points.append(node_radius.get(v, r0) * cp + pos[v])
            case "Bacteria":
                points.append(
                    node_radius.get(v, r0)
                    * np.array([[-1, -1], [1, -1], [0, np.sqrt(3) / 2]])
                    / np.sqrt(2)
                    + pos[v]
                )
            case "Fungi":
                points.append(
                    node_radius.get(v, r0)
                    * np.array([[-1, 1], [1, 1], [-np.sqrt(3) / 2, 0]])
                    / np.sqrt(2)
                    + pos[v]
                )
            case _:
                points.append(node_radius.get(v, r0) * cp + pos[v])

    norm = mc.LogNorm(min(colors), max(colors))

    kwargs.setdefault("facecolors", node_color_map(norm(colors)))

    circles = PolyCollection(points, **inflate_kwargs(H, kwargs))

    ax.add_collection(circles)

    return circles


def _layout_hyper_edges(H, pos, node_radius={}, dr=None):
    """
    Draws a convex hull for each edge in H.

    Position of the nodes in the graph is specified by the position dictionary,
    pos. Convex hulls are spaced out such that if one set contains another, the
    convex hull will surround the contained set. The amount of spacing added
    between hulls is specified by the parameter, dr.

    Parameters
    ----------
    H: Hypergraph
        the entity to be drawn
    pos: dict
        mapping of node and edge positions to R^2
    node_radius: dict
        mapping of node to R^1 (radius of each node)
    dr: float
        the spacing between concentric rings
    ax: Axis
        matplotlib axis on which the plot is rendered

    Returns
    -------
    dict
        A mapping from hyper edge ids to paths (Nx2 numpy matrices)
    """

    if len(node_radius):
        r0 = min(node_radius.values())
    else:
        r0 = get_default_radius(H, pos)

    if dr is None:
        # dr = r0
        dr = 0

    levels = get_set_layering(H)

    radii = {
        v: {v: i for i, v in enumerate(sorted(e, key=levels.get))}
        for v, e in H.dual().edges.elements.items()
    }

    def get_padded_hull(uid, edge):
        # make sure the edge contains at least one node
        if len(edge):
            points = np.vstack(
                [
                    cp * (node_radius.get(v, r0) + dr * (2 + radii[v][uid])) + pos[v]
                    for v in edge
                ]
            )
        # if not, draw an empty edge centered around the location of the edge node (in the bipartite graph)
        else:
            points = 4 * r0 * cp + pos[uid]

        hull = ConvexHull(points)

        return hull.points[hull.vertices]

    return [get_padded_hull(uid, list(H.edges[uid])) for uid in H.edges]


def _draw_hyper_edges(H, pos, ax=None, node_radius={}, dr=None, **kwargs):
    """
    Draws a convex hull around the nodes contained within each edge in H

    Parameters
    ----------
    H: Hypergraph
        the entity to be drawn
    pos: dict
        mapping of node and edge positions to R^2
    node_radius: dict
        mapping of node to R^1 (radius of each node)
    dr: float
        the spacing between concentric rings
    ax: Axis
        matplotlib axis on which the plot is rendered
    kwargs: dict
        keyword arguments, e.g., linewidth, facecolors, are passed through to the PolyCollection constructor

    Returns
    -------
    PolyCollection
        a Matplotlib PolyCollection that can be further styled
    """
    points = _layout_hyper_edges(H, pos, node_radius=node_radius, dr=dr)

    polys = PolyCollection(points, **inflate_kwargs(H.edges, kwargs))

    (ax or plt.gca()).add_collection(polys)

    return polys


##############################
### MICROBIAL HYPERNETWORK ###
##############################
# Back to our custom package code, except the draw method below is also an override of the H.draw method on hypernetx
class Hypernetwork:
    def __init__(
        self,
        population: Population,
        correlation: Correlation,
        hypercorrelation: HyperCorrelation,
        max_group_size: int,
        threshold: float,
    ) -> None:
        self.population = population
        self.correlation = correlation
        self.hypercorrelation = hypercorrelation
        self.max_group_size = max_group_size
        self.threshold = threshold

        self.hyperedge_frames: dict[int, pd.DataFrame] = {}
        self.reduced_hyperedge_frames: dict[int, pd.DataFrame] | None = None
        for n in range(2, max_group_size + 1):
            self.hyperedge_frames[n] = hypercorrelation(
                population=population,
                correlation=correlation,
                group_size=n,
                threshold=threshold,
            )

    def Hypergraph(
        self, reduced: bool = False, minimum_group_size: int = 2
    ) -> Hypergraph:
        """
        Returns a hypernetx hypergraph
        """
        _hyperedge_frames = None
        if reduced:
            if self.reduced_hyperedge_frames is None:
                self.reduce()
            _hyperedge_frames = self.reduced_hyperedge_frames
        else:
            _hyperedge_frames = self.hyperedge_frames

        if _hyperedge_frames is None:
            raise ValueError(
                "_hyperedge_frames is None... this shouldn't be possible. (mostly to make my linter shut up)"
            )

        _hyperedges = []
        for n in range(minimum_group_size, self.max_group_size + 1):
            if _hyperedge_frames.get(n) is None:
                continue

            _hyperedges.extend(_hyperedge_frames[n].index)

        return Hypergraph(_hyperedges)

    def reduce(self) -> None:
        """
        Runs the reduction algorithm on the hyperedges.

        Underneath the class manages two copies of the hyperedges, so this does not modify the existing hyperedge frames.
        """
        if self.reduced_hyperedge_frames is not None:
            LOGGER.debug(
                "A reduce call is being skipped as the reduced_hyperedge_frames is not None."
            )
            return
        try:
            self.reduced_hyperedge_frames = {}
            # The max group size will always be included, so just include those first
            self.reduced_hyperedge_frames[self.max_group_size] = self.hyperedge_frames[
                self.max_group_size
            ].copy()
            for n in range(
                2, self.max_group_size
            ):  # Not + 1, those are always included, see below
                if n not in self.hyperedge_frames:
                    continue

                self.reduced_hyperedge_frames[n] = HyperCorrelation.get_empty_dataframe(
                    n
                )

                for test_hyperedge, test_hypercorrelation in self.hyperedge_frames[
                    n
                ].iterrows():
                    include_flag = True

                    for m in range(n + 1, self.max_group_size + 1):
                        if m not in self.hyperedge_frames:
                            continue

                        for parent_hyperege in self.hyperedge_frames[m].index:
                            if set(test_hyperedge) < set(parent_hyperege):
                                include_flag = False
                                break

                        if not include_flag:
                            break

                    if include_flag:
                        self.reduced_hyperedge_frames[n].loc[
                            test_hyperedge, "Hypercorrelation"
                        ] = test_hypercorrelation["Hypercorrelation"]
        except Exception as e:
            self.reduced_hyperedge_frames = None
            raise e

    def draw(
        self,
        reduced=False,
        minimum_group_size=2,
        edge_color_map=cm["bwr"],
        node_color_map=cm["Greens"],
        pos=None,
        with_color=True,
        with_node_counts=False,
        with_edge_counts=False,
        layout=nx.spring_layout,
        layout_kwargs={},
        dr=None,
        ax=None,
        node_radius=None,
        edges_kwargs={},
        nodes_kwargs={},
        edge_labels_on_edge=True,
        edge_labels={},
        edge_labels_kwargs={},
        node_labels={},
        node_labels_kwargs={},
        with_edge_labels=False,
        with_node_labels=False,
        node_label_alpha=0.35,
        edge_label_alpha=0.35,
        with_additional_edges=None,
        additional_edges_kwargs={},
        return_pos=False,
    ):
        """
        Draw a hypergraph as a Matplotlib figure

        --- CUSTOM DOCS ---
        I'm adding some documentation on arguments or functionality I have modified.

        * reduced is a new bool argument, if true uses the reduced hypergraph.
        * minimum group size can be used to only display larger groups.

        --- HYPERNETX PACKAGE DOCS ---

        By default this will draw a colorful "rubber band" like hypergraph, where
        convex hulls represent edges and are drawn around the nodes they contain.

        This is a convenience function that wraps calls with sensible parameters to
        the following lower-level drawing functions:

        * draw_hyper_edges,
        * draw_hyper_edge_labels,
        * draw_hyper_labels, and
        * draw_hyper_nodes

        The default layout algorithm is nx.spring_layout, but other layouts can be
        passed in. The Hypergraph is converted to a bipartite graph, and the layout
        algorithm is passed the bipartite graph.

        If you have a pre-determined layout, you can pass in a "pos" dictionary.
        This is a dictionary mapping from node id's to x-y coordinates. For example:

            >>> pos = {
            >>> 'A': (0, 0),
            >>> 'B': (1, 2),
            >>> 'C': (5, -3)
            >>> }

        will position the nodes {A, B, C} manually at the locations specified. The
        coordinate system is in Matplotlib "data coordinates", and the figure will
        be centered within the figure.

        By default, this will draw in a new figure, but the axis to render in can be
        specified using :code:`ax`.

        This approach works well for small hypergraphs, and does not guarantee
        a rigorously "correct" drawing. Overlapping of sets in the drawing generally
        implies that the sets intersect, but sometimes sets overlap if there is no
        intersection. It is not possible, in general, to draw a "correct" hypergraph
        this way for an arbitrary hypergraph, in the same way that not all graphs
        have planar drawings.

        Parameters
        ----------
        H: Hypergraph
            the entity to be drawn
        pos: dict
            mapping of node and edge positions to R^2
        with_color: bool
            set to False to disable color cycling of edges
        with_node_counts: bool
            set to True to replace the label for collapsed nodes with the number of elements
        with_edge_counts: bool
            set to True to label collapsed edges with number of elements
        layout: function
            layout algorithm to compute
        layout_kwargs: dict
            keyword arguments passed to layout function
        ax: Axis
            matplotlib axis on which the plot is rendered
        edges_kwargs: dict
            keyword arguments passed to matplotlib.collections.PolyCollection for edges
        node_radius: None, int, float, or dict
            radius of all nodes, or dictionary of node:value; the default (None) calculates radius based on number of collapsed nodes; reasonable values range between 1 and 3
        nodes_kwargs: dict
            keyword arguments passed to matplotlib.collections.PolyCollection for nodes
        edge_labels_on_edge: bool
            whether to draw edge labels on the edge (rubber band) or inside
        edge_labels_kwargs: dict
            keyword arguments passed to matplotlib.annotate for edge labels
        node_labels_kwargs: dict
            keyword argumetns passed to matplotlib.annotate for node labels
        with_edge_labels: bool
            set to False to make edge labels invisible
        with_node_labels: bool
            set to False to make node labels invisible
        node_label_alpha: float
            the transparency (alpha) of the box behind text drawn in the figure for node labels
        edge_label_alpha: float
            the transparency (alpha) of the box behind text drawn in the figure for edge labels
        """

        H = self.Hypergraph(reduced=reduced, minimum_group_size=minimum_group_size)

        ax = ax or plt.gca()

        if pos is None:
            pos = layout_node_link(
                H, with_additional_edges, layout=layout, **layout_kwargs
            )

        r0 = get_default_radius(H, pos)
        a0 = np.pi * r0**2

        def get_node_radius(v):
            if node_radius is None:
                return np.sqrt(a0 * get_collapsed_size(v) / np.pi)
            elif hasattr(node_radius, "get"):
                return node_radius.get(v, 1) * r0
            return node_radius * r0

        # guarantee that node radius is a dictionary mapping nodes to values
        node_radius = {v: get_node_radius(v) for v in H.nodes}

        # for convenience, we are using setdefault to mutate the argument
        # however, we need to copy this to prevent side-effects
        edges_kwargs = edges_kwargs.copy()
        colors = []
        if reduced:
            if self.reduced_hyperedge_frames is None:
                raise ValueError("Somehow no reduced hypergraph exists")
            for edge in getattr(H, "setsystem"):
                colors.append(
                    self.reduced_hyperedge_frames[len(edge)].loc[
                        edge, "Hypercorrelation"
                    ]
                )
        else:
            for edge in getattr(H, "setsystem"):
                colors.append(
                    self.hyperedge_frames[len(edge)].loc[edge, "Hypercorrelation"]
                )

        norm = mc.LogNorm(self.threshold, 1)

        edges_kwargs.setdefault("edgecolors", edge_color_map(norm(colors)))
        edges_kwargs.setdefault("facecolors", "none")

        polys = _draw_hyper_edges(
            H, pos, node_radius=node_radius, ax=ax, dr=dr, **edges_kwargs
        )

        if with_additional_edges:
            nx.draw_networkx_edges(
                with_additional_edges,
                pos=pos,
                ax=ax,
                **inflate_kwargs(
                    with_additional_edges.edges(), additional_edges_kwargs
                ),
            )

        if with_edge_labels:
            labels = get_frozenset_label(
                H.edges, count=with_edge_counts, override=edge_labels
            )

            draw_hyper_edge_labels(
                H,
                pos,
                polys=polys,
                color=edges_kwargs["edgecolors"],
                backgroundcolor=(1, 1, 1, edge_label_alpha),
                labels=labels,
                ax=ax,
                edge_labels_on_edge=edge_labels_on_edge,
                **edge_labels_kwargs,
            )

        if with_node_labels:
            labels = get_frozenset_label(
                H.nodes, count=with_node_counts, override=node_labels
            )

            draw_hyper_labels(
                H,
                pos,
                node_radius=node_radius,
                labels=labels,
                ax=ax,
                va="center",
                xytext=(5, 0),
                textcoords="offset points",
                backgroundcolor=(1, 1, 1, node_label_alpha),
                **node_labels_kwargs,
            )

        _draw_hyper_nodes(
            H,
            pos,
            node_color_map=node_color_map,
            node_radius=node_radius,
            ax=ax,
            **nodes_kwargs,
        )

        if len(H.nodes) == 1:
            x, y = pos[list(H.nodes)[0]]
            s = 20

            ax.axis([x - s, x + s, y - s, y + s])
        else:
            ax.axis("equal")
        ax.axis("off")
        if return_pos:
            return pos


# class _hyperinfo:
#     """
#     Private class to the hypergraph module which helps handle information about hyperedges and their correlation
#     """

#     def __init__(
#         self,
#         hyperedges: dict[str, List[str]] = None,
#         hyperedgecorrs: dict[str, float] = None,
#     ):
#         if hyperedgecorrs is None and hyperedges is None:
#             self.hyperedges = {}
#             self.hyperedgecorrs = {}
#         elif hyperedges is not None and hyperedgecorrs is not None:
#             if set(hyperedges.keys()) != set(hyperedgecorrs.keys()):
#                 raise ValueError(
#                     f"_hyperinfo supplied hyperedges and hyperedgecorrs with different index sets"
#                 )
#             self.hyperedges = hyperedges
#             self.hyperedgecorrs = hyperedgecorrs
#         else:
#             raise ValueError(
#                 "_hyperinfo must be supplied with either nothing or both hyperedges and hyperedgecorrs"
#             )

#     def filter(
#         self, threshold: float = None, max_n: int = None, nodes: List[str] = None
#     ):
#         """
#         Returns a new _hyperinfo instance with filtered information
#         """
#         new_hyperinfo = type(self)()
#         for i in self:
#             if threshold:
#                 if self.hyperedgecorrs[i] < threshold:
#                     continue
#             if max_n:
#                 if len(self.hyperedges[i]) > max_n:
#                     continue
#             if nodes:
#                 remove_flag = True
#                 for node in nodes:
#                     if node in self.hyperedges[i]:
#                         remove_flag = False
#                         break
#                 if remove_flag:
#                     continue

#             # It's passed all filters at this point
#             new_hyperinfo.hyperedges[i] = self.hyperedges[i].copy()
#             new_hyperinfo.hyperedgecorrs[i] = self.hyperedgecorrs[i]

#         return new_hyperinfo

#     def __getitem__(self, key):
#         """
#         Convenience method to get both the edge and the correlation
#         """
#         return (self.hyperedges[key], self.hyperedgecorrs[key])

#     def pop(self, key):
#         """
#         pops the hyperedge and hyperedgecorr associated to the index
#         """
#         return (self.hyperedges.pop(key), self.hyperedgecorrs.pop(key))

#     def __iter__(self):
#         return iter(self.hyperedges)


# def edgecount(self):
#     """
#     Returns a pandas dataframe whose index are otus and whose values are the number of hyperedges that that otu is contained in
#     """
#     _d = {}
#     for i in self._hyperinfo:
#         for otu in self._hyperinfo.hyperedges[i]:
#             if otu in _d:
#                 _d[otu] += 1
#             else:
#                 _d[otu] = 1
#     return pd.DataFrame(
#         list(_d.values()), index=list(_d.keys()), columns=["Hyperedge Count"]

# def export(
#     self,
#     export_dir: str,
#     cmap=cm["copper"],
#     with_edge_labels=False,
#     with_node_labels=False,
# ) -> None:
#     """
#     Exports a variety of information about this particular hypergraph.
#     """

#     export_dir = self.get_export_dir(export_dir)

#     if len(self._hyperinfo.hyperedges) == 0:
#         return

#     H = hnx.Hypergraph(self._hyperinfo.hyperedges)
#     H_stats = hnx_stats.dist_stats(H)

#     ### Summary tab ###
#     summary = pd.DataFrame(
#         columns=["Value"], index=["Total OTUs", "Total Hyperedges"]
#     )
#     summary.loc["Total Hyperedges", "Value"] = H_stats["ncols"]
#     summary.loc["Total OTUs", "Value"] = H_stats["nrows"]
#     summary.loc["Density", "Value"] = H_stats["density"]
#     hypercorrs = list(self._hyperinfo.hyperedgecorrs.values())
#     summary.loc["Max Hypercorrelation", "Value"] = np.max(hypercorrs)
#     summary.loc["Min Hypercorrelation", "Value"] = np.min(hypercorrs)
#     summary.loc["Average Hypercorrelation", "Value"] = np.mean(hypercorrs)
#     summary.loc["Median Hypercorrelation", "Value"] = np.median(hypercorrs)
#     summary.loc["Component Sizes", "Value"] = str(hnx_stats.comp_dist(H))

#     ### Write to the summary excel
#     with pd.ExcelWriter(os.path.join(export_dir, "data.xlsx")) as writer:
#         summary.to_excel(writer, sheet_name="Summary")
#         self.edgecount().sort_values("Hyperedge Count", ascending=False).to_excel(
#             writer, sheet_name="Edge Counts"
#         )
#         for i, df in self.data_frames.items():
#             df.to_excel(writer, sheet_name=f"{i}-data")

#     l = list(self._hyperinfo.hyperedgecorrs.values())
#     norm = plt.Normalize(min(l), max(l))
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm.set_array([])
#     colors = sm.to_rgba(l)

#     self.draw(
#         H,
#         with_edge_labels=with_edge_labels,
#         with_node_labels=with_node_labels,
#         dr=10**-3,
#         edges_kwargs={"edgecolors": colors},
#         nodes_kwargs={},
#     )
#     plt.savefig(os.path.join(export_dir, "drawing.png"))
