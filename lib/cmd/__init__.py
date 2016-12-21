import importlib

commands = ['gengrid', 'preproc', 'maplatlon']


def get_module(name):
    return importlib.import_module("%s.%s" % (__name__, name))
