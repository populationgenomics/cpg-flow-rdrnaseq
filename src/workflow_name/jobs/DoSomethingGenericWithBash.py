"""
A job should contain the logic for a single Stage
"""

from cpg_utils.hail_batch import get_batch
from cpg_utils.config import config_retrieve
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def echo_statement_to_file(statement: str, output_file: str) -> 'Job':
    """
    This is a simple example of a job that writes a statement to a file.

    Args:
        statement (str): the intended file contents
        output_file (str): the path to write the file to

    Returns:
        the resulting job
    """

    # create a job
    j = get_batch().new_job(f'echo "{statement}" to {output_file}')

    # choose an image to run this job in (default is bare ubuntu)
    j.image(config_retrieve(['workflow', 'driver_image']))

    # write the statement to the file
    j.command(f'echo "{statement}" > {j.output}')

    # write the output to the expected location
    get_batch().write_output(j.output, output_file)

    # return the job
    return j
