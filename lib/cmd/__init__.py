import importlib

commands = ['gengrid', 'preproc', 'projlatlon', 'calcweights']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
