import importlib

commands = ['gengrid', 'preproc', 'projlatlon', 'calcweights', 'applyweights']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
