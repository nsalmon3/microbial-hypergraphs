"""
Manages hypergraphs
"""

from hypernetx import Hypergraph


class Hypernetwork:
    def __init__(self, hyperedges: dict[str, list[str]]) -> None:
        self.hyperedges = hyperedges

    @property
    def H(self) -> Hypergraph:
        return Hypergraph(self.hyperedges)
