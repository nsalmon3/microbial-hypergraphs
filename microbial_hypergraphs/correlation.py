"""
Manages simple correlations
"""

from typing import Any, Callable

import pandas as pd

from microbial_hypergraphs.core.exception import MicrobialHypergraphException
from microbial_hypergraphs.core import LOGGER
from microbial_hypergraphs.population import Population


##################
### EXCEPTIONS ###
##################
class CorrelationNameNotFoundException(MicrobialHypergraphException):
    def __init__(self, not_found_name: str) -> None:
        super().__init__("Could not find a correlation with name %s", not_found_name)


############
### TYPE ###
############
class Correlation:
    """
    Manages a method of correlating sample data with itself
    """

    def __init__(
        self,
        name: str,
        description: str,
        correlator: Callable[[Population], pd.DataFrame],
    ):
        self._correlator = correlator
        self.name = name
        self.description = description

        LOGGER.debug("Initialized correlation %s.", self)

    def __call__(self, population: Population) -> Any:
        return self._correlator(population)

    def __str__(self) -> str:
        return self.name


#################
### INSTANCES ###
#################
CORRELATION_INSTANCES: dict[str, Correlation] = {}


# Pearson Correlation
def _pearson_correlator(population: Population) -> pd.DataFrame:
    return population.otu_samples.corr("pearson")


CORRELATION_INSTANCES["Pearson"] = Correlation(
    "Pearson", "The standard Pearson correlation.", _pearson_correlator
)
