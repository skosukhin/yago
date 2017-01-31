import importlib

commands = ['gengrid', 'project', 'overlap', 'append', 'calcmask',
            'applymask', 'slice', 'interp']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
