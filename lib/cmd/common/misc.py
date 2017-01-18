import os


def create_dir_for_file(filename):
    try:
        os.makedirs(os.path.dirname(filename))
    except:
        pass



