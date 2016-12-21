import importlib

commands = ['gengrid', 'preproc', 'projlatlon']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
