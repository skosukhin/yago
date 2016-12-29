import importlib

commands = ['gengrid', 'preproc', 'projlatlon', 'calcweights', 'applyweights',
            'overlap', 'append']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
