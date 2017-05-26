import numpy as np


class Translator(object):
    """
    Applies false easting and false northing.
    """

    def __init__(self, easting, northing):
        self._easting = easting
        self._northing = northing

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self._easting == other._easting and
                    self._northing == other._northing)
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self._easting, self._northing))

    @property
    def easting(self):
        return self._easting

    @property
    def northing(self):
        return self._northing

    def convert_points(self, xx, yy):
        return np.asanyarray(xx) + self._easting, \
               np.asanyarray(yy) + self._northing

    def restore_points(self, xx, yy):
        return np.asanyarray(xx) - self._easting, \
               np.asanyarray(yy) - self._northing
