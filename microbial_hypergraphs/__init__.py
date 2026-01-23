"""
Initializes the microbial_hypergraphs pacakge
"""

from pathlib import Path

from microbial_hypergraphs.core import DATA_PATH

from microbial_hypergraphs.population import POPULATION_INSTANCES
from microbial_hypergraphs.correlation import CORRELATION_INSTANCES

from microbial_hypergraphs.api import *

############
### DATA ###
############
for population_instance in POPULATION_INSTANCES:
    if not Path.exists(DATA_PATH / "cache" / population_instance):
        Path.mkdir(DATA_PATH / "cache" / population_instance, parents=True)

    for correlation_instance in CORRELATION_INSTANCES:
        if not Path.exists(
            DATA_PATH / "cache" / population_instance / correlation_instance
        ):
            Path.mkdir(DATA_PATH / "cache" / population_instance / correlation_instance)
