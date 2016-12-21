
class Projector(object):
    short_name = None
    long_name = None
    standard_name = None

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def convert_point(self, la, lo):
        raise NotImplementedError()

    def restore_point(self, x, y):
        raise NotImplementedError()

    def convert_vector(self, u, v, la, lo, return_point=False):
        raise NotImplementedError()

    def restore_vector(self, u, v, x, y, return_point=False):
        raise NotImplementedError()
