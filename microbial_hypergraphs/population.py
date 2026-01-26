"""
Defines the Population class, used to represent different populations of OTUs
"""

import pandas as pd

from microbial_hypergraphs.core.data import get_otu_samples, get_sample_info
from microbial_hypergraphs.core.exception import MicrobialHypergraphException
from microbial_hypergraphs.core import LOGGER

_all_otu_samples = get_otu_samples()
_sample_info = get_sample_info()


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

    def __init__(
        self, name: str, description: str, otus: list[str], samples: list[int]
    ):
        self._validate_otu_list(otus)
        self._validate_sample_list
        self.otus = otus
        self.samples = samples
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

    @classmethod
    def _validate_sample_list(cls, sample_list: list[int]) -> None:
        for sample in sample_list:
            if sample not in _all_otu_samples.index:
                raise Exception(f"{sample} not in list of samples.")

    @property
    def otu_samples(self) -> pd.DataFrame:
        """
        A dataframe representing this populations sample values

        :return: DataFrame representing this population
        :rtype: DataFrame
        """
        return _all_otu_samples.loc[self.samples, self.otus]


#################
### INSTANCES ###
#################
POPULATION_INSTANCES: dict[str, Population] = {}

POPULATION_INSTANCES["PresentInAllSamples"] = Population(
    name="PresentInAllSamples",
    description="OTUs that show up in ALL samples.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
    samples=_sample_info.index,
)

# Present in all samples for each country
POPULATION_INSTANCES["PresentInAllSamples_Spain"] = Population(
    name="PresentInAllSamples_Spain",
    description="OTUs that show up in ALL samples. Only Spain samples considered for statistics.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
    samples=_sample_info.loc[_sample_info["Country"] == "Spain"].index,
)
POPULATION_INSTANCES["PresentInAllSamples_France"] = Population(
    name="PresentInAllSamples_France",
    description="OTUs that show up in ALL samples. Only France samples considered for statistics.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
    samples=_sample_info.loc[_sample_info["Country"] == "France"].index,
)
POPULATION_INSTANCES["PresentInAllSamples_Sweden"] = Population(
    name="PresentInAllSamples_Sweden",
    description="OTUs that show up in ALL samples. Only Sweden samples considered for statistics.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
    samples=_sample_info.loc[_sample_info["Country"] == "Sweden"].index,
)
POPULATION_INSTANCES["PresentInAllSamples_Germany"] = Population(
    name="PresentInAllSamples_Germany",
    description="OTUs that show up in ALL samples. Only Germany samples considered for statistics.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
    samples=_sample_info.loc[_sample_info["Country"] == "Germany"].index,
)
POPULATION_INSTANCES["PresentInAllSamples_Switzerland"] = Population(
    name="PresentInAllSamples_Sweden",
    description="OTUs that show up in ALL samples. Only Switzerland samples considered for statistics.",
    otus=_all_otu_samples.loc[:, (_all_otu_samples > 0).all(axis=0)].columns.tolist(),
    samples=_sample_info.loc[_sample_info["Country"] == "Switzerland"].index,
)
