class Converter(object):
    """
    Links the geographical coordinate system with a Cartesian coordinate
    system by combining rotation and projection transformations.
    """

    def __init__(self, rotor, projection, translator):
        """
        Constructor of the class.
        :param rotor: Instance of the class Rotor
        :param projection: Instance of the class Projection.
        """
        self._rotor = rotor
        self._projection = projection
        self._translator = translator

    @property
    def rotor(self):
        return self._rotor

    @property
    def projection(self):
        return self._projection

    @property
    def translator(self):
        return self._translator

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def convert_points(self, lats, lons):
        rot_lats, rot_lons = self._rotor.convert_points(lats, lons)
        xx, yy = self._projection.convert_points(rot_lats, rot_lons)
        xx, yy = self._translator.convert_points(xx, yy)
        return xx, yy

    def restore_points(self, xx, yy):
        xx, yy = self._translator.restore_points(xx, yy)
        rot_lats, rot_lons = self._projection.restore_points(xx, yy)
        lats, lons = self._rotor.restore_points(rot_lats, rot_lons)
        return lats, lons

    def convert_vectors(self, uu, vv, lats, lons, return_points=False):
        rot_uu, rot_vv, rot_lats, rot_lons = \
            self._rotor.convert_vectors(uu, vv, lats, lons, True)
        if not return_points:
            return self._projection.convert_vectors(rot_uu, rot_vv,
                                                    rot_lats, rot_lons,
                                                    False)
        else:
            proj_uu, proj_vv, proj_xx, proj_yy = \
                self._projection.convert_vectors(rot_uu, rot_vv, rot_lats,
                                                 rot_lons, True)
            xx, yy = self._translator.convert_points(proj_xx, proj_yy)
            return proj_uu, proj_vv, xx, yy

    def restore_vectors(self, uu, vv, xx, yy, return_points=False):
        if return_points:
            xx, yy = self._translator.restore_points(xx, yy)
        rot_uu, rot_vv, rot_lats, rot_lons = \
            self._projection.restore_vectors(uu, vv, xx, yy, True)
        return self._rotor.restore_vectors(rot_uu, rot_vv, rot_lats, rot_lons,
                                           return_points)
