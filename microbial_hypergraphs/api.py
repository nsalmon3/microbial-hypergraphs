"""
Exposed objects when importing the package
"""

import pandas as pd

import rich.table

from microbial_hypergraphs.correlation import (
    CORRELATION_INSTANCES,
    CorrelationNameNotFoundException,
)
from microbial_hypergraphs.population import (
    POPULATION_INSTANCES,
    PopulationNameNotFoundException,
)
from microbial_hypergraphs.hypercorrelation import (
    HYPERCORRELATION_INSTANCES,
    HypercorrelationNameNotFoundException,
)
from microbial_hypergraphs.core import CONSOLE


def hypercorrelate(
    population_name: str,
    correlation_name: str,
    hypercorrelation_name: str,
    group_size: int,
    threshold: float,
) -> pd.DataFrame:
    if hypercorrelation_name not in HYPERCORRELATION_INSTANCES:
        raise HypercorrelationNameNotFoundException(hypercorrelation_name)
    if population_name not in POPULATION_INSTANCES:
        raise PopulationNameNotFoundException(population_name)
    if correlation_name not in CORRELATION_INSTANCES:
        raise CorrelationNameNotFoundException(correlation_name)

    return HYPERCORRELATION_INSTANCES[hypercorrelation_name](
        threshold=threshold,
        group_size=group_size,
        population=POPULATION_INSTANCES[population_name],
        correlation=CORRELATION_INSTANCES[correlation_name],
    )


def print_names():
    populations = rich.table.Table(title="Populations")
    populations.add_column("Name")
    populations.add_column("Description")

    for population in POPULATION_INSTANCES.values():
        populations.add_row(population.name, population.description)

    correlations = rich.table.Table(title="Correlations")
    correlations.add_column("Name")
    correlations.add_column("Description")

    for correlation in CORRELATION_INSTANCES.values():
        correlations.add_row(correlation.name, correlation.description)

    hypercorrelations = rich.table.Table(title="Hypercorrelations")
    hypercorrelations.add_column("Name")
    hypercorrelations.add_column("Description")

    for hypercorrelation in HYPERCORRELATION_INSTANCES.values():
        hypercorrelations.add_row(hypercorrelation.name, hypercorrelation.description)

    CONSOLE.print(populations)
    CONSOLE.print(correlations)
    CONSOLE.print(hypercorrelations)
