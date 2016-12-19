import importlib

commands = ['gengrid', 'strip', 'maplatlon']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
