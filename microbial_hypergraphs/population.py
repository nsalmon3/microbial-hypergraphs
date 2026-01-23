"""
Defines the Population class, used to represent different populations of OTUs
"""

import pandas as pd

from microbial_hypergraphs.core.data import get_otu_samples
from microbial_hypergraphs.core.exception import MicrobialHypergraphException
from microbial_hypergraphs.core import LOGGER

_all_otu_samples = get_otu_samples()


##################
### EXCEPTIONS ###
##################
class PopulationNameNotFoundException(MicrobialHypergraphException):
    def __init__(self, not_found_name: str) -> None:
        super().__init__("Could not find a population with name %s", not_found_name)


############
### TYPE ###
############
class Population:
    """
    Represents a particular population of OTUs
    """

    def __init__(self, name: str, description: str, otus: list[str]):
        self._validate_otu_list(otus)
        self.otus = otus
        self.name = name
        self.description = description

        LOGGER.debug("Initalized population %s.", self)

    def __str__(self) -> str:
        return self.name

    @classmethod
    def _validate_otu_list(cls, otu_list: list[str]) -> None:
        """
        Checks that the list of strings only contains valid OTU_ids
        """
        for otu in otu_list:
            if otu not in _all_otu_samples.columns:
                raise Exception(f"{otu} not in list of OTUs")

    @property
    def otu_samples(self) -> pd.DataFrame:
        """
        A dataframe representing this populations sample values

        :return: DataFrame representing this population
        :rtype: DataFrame
        """
        return _all_otu_samples[self.otus]


#################
### INSTANCES ###
#################
POPULATION_INSTANCES: dict[str, Population] = {}
POPULATION_INSTANCES["PresentInAllSamples"] = Population(
    name="PresentInAllSamples",
    description="OTUs that show up in ALL samples.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
)
