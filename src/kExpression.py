import sys
import click

from src.click_context import cli
from src.batch_query import batchQuery      # pylint: disable=relative-beyond-top-level

# cli.add_command(index_main, name="index")
cli.add_command(batchQuery, name="batchQuery")

if __name__ == '__main__':
    cli()