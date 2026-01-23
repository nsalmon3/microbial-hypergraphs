"""
Command Line Interface
"""

from typing import Annotated

import typer

from microbial_hypergraphs.api import hypercorrelate, print_names
from microbial_hypergraphs.core import LOGGER, EXPORT_PATH

cli = typer.Typer()


@cli.command(name="list")
def cli_list():
    """
    Lists the relevant objects for cli arguments
    """
    print_names()


@cli.command(name="hypercorrelate")
def cli_hypercorrelate(
    population_name: Annotated[str, typer.Option(..., help="Name of the Population.")],
    correlation_name: Annotated[
        str, typer.Option(..., help="Name of the correlation method.")
    ],
    hypercorrelation_name: Annotated[
        str, typer.Option(..., help="Name of the hypercorrelation method.")
    ],
    group_size: Annotated[int, typer.Option(..., help="Number of OTUs in each group.")],
    threshold: Annotated[
        float,
        typer.Option(
            ...,
            help="Minimal hypercorrelation value for a group to be included in the result.",
        ),
    ],
):
    """
    Creates a hypercorrelation dataframe and exports it as a csv
    """

    df = hypercorrelate(
        population_name=population_name,
        correlation_name=correlation_name,
        hypercorrelation_name=hypercorrelation_name,
        group_size=group_size,
        threshold=threshold,
    )

    export_filename = (
        EXPORT_PATH
        / f"{population_name}_{correlation_name}_{hypercorrelation_name}_{int(100 * threshold)}_{group_size}.csv"
    )
    df.to_csv(export_filename)

    LOGGER.info("Exported results to %s", export_filename)
