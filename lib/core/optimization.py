import numpy as np


def great_circle_distances(lats_deg0, lons_def0,
                           lats_deg1, lons_deg1,
                           earth_radius):
    lats_rad0 = np.radians(lats_deg0)
    lats_rad1 = np.radians(lats_deg1)
    half_delta_lons_rad = np.radians(lons_deg1 - lons_def0) / 2.0
    half_delta_lats_rad = (lats_rad1 - lats_rad0) / 2.0

    a = (np.power(np.sin(half_delta_lats_rad), 2.0) +

         np.cos(lats_rad0) * np.cos(lats_rad1) *
         np.power(np.sin(half_delta_lons_rad), 2.0))

    return earth_radius * 2.0 * np.arcsin(np.sqrt(a))


def grid_distances(lats_deg, lons_deg, earth_radius):
    along_x = great_circle_distances(lats_deg[:-1], lons_deg[:-1],
                                     lats_deg[1:], lons_deg[1:],
                                     earth_radius)
    along_y = great_circle_distances(lats_deg[:, :-1], lons_deg[:, :-1],
                                     lats_deg[:, 1:], lons_deg[:, 1:],
                                     earth_radius)
    return along_x, along_y

