"""
mini script which prints "Hello, world!"
if an argument is provided, instead prints "Hello, {arg}!"
"""

from argparse import ArgumentParser


def main(name_string: str):
    """
    Main function, which prints a greeting
    """
    print(f'Hello, {name_string}!')


def cli_main():
    """
    CLI entrypoint, if this is run as a script
    """
    parser = ArgumentParser()

    # take one optional positional argument
    parser.add_argument(
        'name',
        help='Name to greet',
        default='world',
        nargs='?',  # optional positional argument
        type=str,
    )
    args = parser.parse_args()

    main(args.name)


if __name__ == '__main__':
    cli_main()
