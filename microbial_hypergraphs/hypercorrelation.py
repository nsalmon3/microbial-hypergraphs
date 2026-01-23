"""
Manages hypercorrelations
"""

from itertools import combinations
from typing import Callable, Iterator
import time
from datetime import timedelta, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import math

from microbial_hypergraphs.core.exception import MicrobialHypergraphException
from microbial_hypergraphs.population import Population
from microbial_hypergraphs.correlation import Correlation
from microbial_hypergraphs.core import LOGGER, DATA_PATH


##################
### EXCEPTIONS ###
##################
class HypercorrelationNameNotFoundException(MicrobialHypergraphException):
    def __init__(self, not_found_name: str) -> None:
        super().__init__(
            "Could not find a hyeprcorrelation with name %s", not_found_name
        )


############
### TYPE ###
############
class HyperCorrelation:
    """
    Manages a specific hypercorrelation
    """

    def __init__(
        self,
        name: str,
        description: str,
        hypercorrelator: Callable[[float, int, Population, Correlation], pd.DataFrame],
    ) -> None:
        self._hypercorrelator = hypercorrelator
        self.name = name
        self.description = description

        LOGGER.debug("Initalized hypercorrelation %s.", self)

    def __call__(
        self,
        population: Population,
        correlation: Correlation,
        group_size: int,
        threshold: float,
    ) -> pd.DataFrame:
        cache_file_path = self.cache_path(
            population=population, correlation=correlation
        )
        # First we check whether this information has already been computed
        # The information on what computations for a particular hypercorrelation will be stored in the log sheet
        logs = None
        if Path.exists(cache_file_path):
            LOGGER.info("Checking cached hypercorrelation results.")
            logs = pd.read_excel(cache_file_path, sheet_name="logs", index_col=0)

            if len(logs.loc[logs["group_size"] == group_size]) > 0:
                min_threshold_for_this_group_size = min(
                    logs.loc[logs["group_size"] == group_size]["threshold"]
                )
            else:
                min_threshold_for_this_group_size = np.inf

            if min_threshold_for_this_group_size <= threshold:
                # This means that we must have already computed for this particular hypercorrelation
                LOGGER.info(
                    "Cached results already contain equal or lower hypercorrelations (%f) than threshold, using cached result.",
                    min_threshold_for_this_group_size,
                )
                _df = pd.read_excel(
                    cache_file_path,
                    sheet_name=str(group_size),
                    index_col=list(range(group_size)),
                )
                return _df.loc[_df["Hypercorrelation"] >= threshold].copy()
            else:
                LOGGER.info(
                    "Minimal hypercorrelation found in cached results was %f, but passed threshold is %f, will not use cached hypercorrelation.",
                    min_threshold_for_this_group_size,
                    threshold,
                )

        LOGGER.info(
            "Beginning a %s hypercorrelation with parameters:\n\tthreshold: %f\n\tgroup_size: %i\n\tpopulation: %s\n\tcorrelation: %s",
            self,
            threshold,
            group_size,
            population,
            correlation,
        )
        start_time = time.perf_counter()
        return_value = self._hypercorrelator(
            threshold, group_size, population, correlation
        )
        LOGGER.info(
            "Successfully found %i out of %i groups were above threshold %f in %s.",
            len(return_value.index),
            math.comb(len(population.otus), group_size),
            threshold,
            str(timedelta(seconds=time.perf_counter() - start_time)),
        )

        # Now cache the result
        with pd.ExcelWriter(
            cache_file_path,
            mode="a" if Path.exists(cache_file_path) else "w",
            if_sheet_exists="replace" if Path.exists(cache_file_path) else None,
        ) as writer:
            return_value.to_excel(writer, sheet_name=str(group_size))

        # Next log that it happened
        new_log = pd.DataFrame(
            index=[datetime.now()],
            columns=["group_size", "threshold"],
            data=[[group_size, threshold]],
        )
        if logs is not None:
            with pd.ExcelWriter(
                cache_file_path, mode="a", if_sheet_exists="replace"
            ) as writer:
                pd.concat(
                    [
                        logs,
                        new_log,
                    ]
                ).to_excel(
                    writer,
                    sheet_name="logs",
                )
        else:
            with pd.ExcelWriter(cache_file_path, mode="a") as writer:
                new_log.to_excel(writer, sheet_name="logs")

        LOGGER.info("Successfully cached result to %s", cache_file_path)

        return return_value

    def __str__(self) -> str:
        return self.name

    @classmethod
    def get_iter(cls, group_size: int, population: Population) -> Iterator:
        return combinations(population.otus, group_size)

    @classmethod
    def get_empty_dataframe(cls, group_size: int) -> pd.DataFrame:
        return pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [[] for _ in range(1, group_size + 1)],
                names=[f"OTU {n}" for n in range(1, group_size + 1)],
            ),
            columns=["Hypercorrelation"],
            dtype=float,
        )

    def cache_path(self, population: Population, correlation: Correlation) -> Path:
        return Path(
            DATA_PATH / "cache" / population.name / correlation.name / f"{self}.xlsx"
        )


#################
### INSTANCES ###
#################

HYPERCORRELATION_INSTANCES: dict[str, HyperCorrelation] = {}


# Minimal Hypercorrelation
#   The Minimal Hypercorrelation takes the smallest value amongst all correlation values
def _minimal_hypercorrelator(
    threshold: float, group_size: int, population: Population, correlation: Correlation
) -> pd.DataFrame:
    _df = HyperCorrelation.get_empty_dataframe(group_size=group_size)
    _iterator = HyperCorrelation.get_iter(group_size=group_size, population=population)

    correlation_df = correlation(population=population)

    for group in _iterator:
        _min = np.inf
        for pair in combinations(group, 2):
            _min = min(_min, correlation_df.loc[pair[0], pair[1]])
        if _min >= threshold:
            _df.loc[group] = _min

    return _df


# Arithmetic Hypercorrelation
#   The Arithmetic Hypercorrelation takes the average value amongst all correlation values
def _arithmetic_hypercorrelator(
    threshold: float, group_size: int, population: Population, correlation: Correlation
) -> pd.DataFrame:
    _df = HyperCorrelation.get_empty_dataframe(group_size=group_size)
    _iterator = HyperCorrelation.get_iter(group_size=group_size, population=population)

    correlation_df = correlation(population=population)

    _group_size = math.comb(group_size, 2)
    for group in _iterator:
        _sum = 0
        for pair in combinations(group, 2):
            _sum += correlation_df.loc[pair[0], pair[1]]
        if _sum / _group_size >= threshold:
            _df.loc[group] = _sum / _group_size

    return _df


HYPERCORRELATION_INSTANCES["Minimal"] = HyperCorrelation(
    "Minimal",
    "The smallest pairwise correlation in the group.",
    _minimal_hypercorrelator,
)
HYPERCORRELATION_INSTANCES["Arithmetic"] = HyperCorrelation(
    "Arithmetic",
    "The arithmetic mean of the pairwise correlations in the group.",
    _arithmetic_hypercorrelator,
)
