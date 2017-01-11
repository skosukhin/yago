from core.projections.lambert import LambertConformalProjection
from core.projections.mercator import MercatorProjection
from core.projections.polar_stereographic import PolarStereographicProjection
from core.projections.sinusoidal import SinusoidalProjection

_proj_list = [LambertConformalProjection, MercatorProjection,
              PolarStereographicProjection, SinusoidalProjection]

projections = dict(
    [(proj_cls.short_name, proj_cls) for proj_cls in _proj_list])
