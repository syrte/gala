# coding: utf-8

""" Astropy coordinate class for the Sagittarius coordinate system """

from __future__ import division, print_function


# Third-party
import numpy as np

import astropy.coordinates as coord
from astropy.coordinates import frame_transform_graph
try:
    from astropy.coordinates.matrix_utilities import matrix_transpose
    ASTROPY_1_3 = True
except ImportError:
    from .matrix_utilities import matrix_transpose
    ASTROPY_1_3 = False

if not ASTROPY_1_3:
    import astropy
    import warnings
    warnings.warn("We recommend using Astropy v1.3 or later. You have: {}"
                  .format(astropy.__version__), DeprecationWarning)

__all__ = ["GD1"]

class GD1(coord.BaseCoordinateFrame):
    """
    A Heliocentric spherical coordinate system defined by the orbit
    of the GD1 stream, as described in
    Koposov et al. 2010 (see: `<http://arxiv.org/abs/0907.1085>`_).

    For more information about this class, see the Astropy documentation
    on coordinate frames in :mod:`~astropy.coordinates`.

    Parameters
    ----------
    representation : :class:`~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    phi1 : angle_like, optional, must be keyword
        The longitude-like angle corresponding to Orphan's orbit.
    phi2 : angle_like, optional, must be keyword
        The latitude-like angle corresponding to Orphan's orbit.
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    """
    default_representation = coord.SphericalRepresentation

    frame_specific_representation_info = {
        'spherical': [coord.RepresentationMapping('lon', 'phi1'),
                      coord.RepresentationMapping('lat', 'phi2'),
                      coord.RepresentationMapping('distance', 'distance')],
        'unitspherical': [coord.RepresentationMapping('lon', 'phi1'),
                          coord.RepresentationMapping('lat', 'phi2')]
    }

# Rotation matrix as defined in the Appendix of Koposov et al. (2010)
R = np.array([[-0.4776303088, -0.1738432154, 0.8611897727],
              [0.510844589, -0.8524449229, 0.111245042],
              [0.7147776536, 0.4930681392, 0.4959603976]])

@frame_transform_graph.transform(coord.StaticMatrixTransform, coord.ICRS, GD1)
def icrs_to_gd1():
    """ Compute the transformation from Galactic spherical to
        heliocentric GD1 coordinates.
    """
    return R

@frame_transform_graph.transform(coord.StaticMatrixTransform, GD1, coord.ICRS)
def gd1_to_icrs():
    """ Compute the transformation from heliocentric GD1 coordinates to
        spherical Galactic.
    """
    return matrix_transpose(icrs_to_gd1())
