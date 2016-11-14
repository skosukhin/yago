import importlib

commands = ['gengrid']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
