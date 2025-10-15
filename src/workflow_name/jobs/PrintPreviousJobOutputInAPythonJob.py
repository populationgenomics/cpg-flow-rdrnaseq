"""
trivially simple python job example, using a utility constant
"""

from workflow_name.utils import DATE_STRING


def print_file_contents(input_file: str) -> str:
    """
    This is a simple example of a job that prints the contents of a file.

    Args:
        input_file (str): the path to the file to print
    """

    with open(input_file) as f:
        contents = f.read()

    print(f'Contents of {input_file} on {DATE_STRING}:')
    print(contents)
    return contents
