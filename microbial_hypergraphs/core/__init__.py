"""
Initializes core functionality
"""

from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

from rich.console import Console
from rich.logging import RichHandler

from microbial_hypergraphs.core.env import ROOT_PATH
from microbial_hypergraphs.core.data import DATA_PATH

START_TIMESTAMP = datetime.now()

###############
### EXPORTS ###
###############
EXPORT_PATH = ROOT_PATH / "exports" / str(START_TIMESTAMP.year) / str(START_TIMESTAMP.month) / str(START_TIMESTAMP.day)
Path.mkdir(EXPORT_PATH, parents=True) if not Path.exists(EXPORT_PATH) else None

###############
### LOGGING ###
###############

# Create logger
LOGGER = logging.getLogger('core')
LOGGER.setLevel(logging.DEBUG)

# File logging setup
_logging_dir = ROOT_PATH / "logs" / str(START_TIMESTAMP.year) / str(START_TIMESTAMP.month) / str(START_TIMESTAMP.day)
Path.mkdir(_logging_dir, parents=True) if not Path.exists(_logging_dir) else None

_file_handler = RotatingFileHandler(
    filename= _logging_dir / START_TIMESTAMP.strftime("%H_%M_%S.log")
)
_file_handler.setLevel(logging.DEBUG)
_file_handler.setFormatter(
    fmt=logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(message)s"
    )
)

LOGGER.addHandler(hdlr=_file_handler)

###############
### CONSOLE ###
###############

CONSOLE = Console()

_console_handler = RichHandler(
    level=logging.INFO,
    console=CONSOLE
)
_console_handler.setFormatter(fmt=logging.Formatter(fmt="%(message)s"))
LOGGER.addHandler(hdlr=_console_handler)

############
### DATA ###
############
if not Path.exists(DATA_PATH / "hypercorrelation"):
    Path.mkdir(DATA_PATH / "hypercorrelation", parents=True)

