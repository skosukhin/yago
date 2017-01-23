import importlib

description = 'performs bilinear interpolation of fields specified on a ' \
              'structured grid in Cartesian coordinates'

interp_commands = ['calc', 'apply']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))


def setup_parser(parser):
    subparsers = parser.add_subparsers(metavar='command',
                                       dest='interp_command')
    for c in interp_commands:
        module = get_module(c)
        sub = subparsers.add_parser(c, help=module.description)
        module.setup_parser(sub)


def cmd(args):
    get_module(args.interp_command).cmd(args)
