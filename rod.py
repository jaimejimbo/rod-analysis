"""
Rod object.
It has needed methods and variables.
"""

import math
from methods import is_in_circle



class Rod(object):
    """
    Rod object.
    """

    def __init__(self, (ID, area, xm, ym, major, minor,
                        angle, feret, feretx, ferety,
                        feretangle, minferet, xstart, ystart)):
        """
        Initialization of rod
        """                                     #Column
        self._id = int(ID)                      #0
        self._area = float(area)                #1
        self._x_mid = float(xm)                 #2
        self._y_mid = float(ym)                 #3
        self._major = float(major)              #4
        self._minor = float(minor)              #5
        self._angle = float(angle)              #6
        self._feret = float(feret)              #7
        self._feret_x = float(feretx)           #8
        self._feret_y = float(ferety)           #9
        self._feret_angle = float(feretangle)   #10
        self._min_feret = float(minferet)       #11
        self._x_start = float(xstart)           #12
        self._y_start = float(ystart)           #13
        self._hash = 0
        self._kappa = float(self.feret)/self.min_feret

    @property
    def feret(self):
        """
        Feret length.
        (wikipedia) The Feret diameter or Feret's diameter is a measure
        of an object size along a specified direction. In general, it can
        be defined as the distance between the two parallel planes restricting
        the object perpendicular to that direction. It is therefore also called
        the caliper diameter, referring to the measurement of the object size
        with a caliper. This measure is used in the analysis of particle sizes,
        for example in microscopy, where it is applied to projections of a
        three-dimensional (3D) object on a 2D plane. In such cases, the Feret
        diameter is defined as the distance between two parallel tangential
        lines rather than planes.[1][2]
        """
        return self._feret

    def __eq__(self, rod2):
        """
        Check if a rod is the same as another rod.
        Rods must be of the same group.
        """
        return self.hash_ == rod2.hash_

    def __ne__(self, rod2):
        """
        != magic method
        """
        return not self == rod2

    @property
    def identifier(self):
        """
        Returns an identification number.
        """
        return self._id

    @property
    def hash_(self):
        """
        Returns an unique number of this rod.
        Uses some of rod properties.
        """
        output = ""
        output += str(self.identifier)
        output += str(int(self.min_feret))
        output += str(int(self.x_mid))
        output += str(int(self.y_mid))
        output += str(int(self.kappa))
        return int(output)


    @property
    def min_feret(self):
        """
        Minimum Feret length.
        """
        return self._min_feret

    @property
    def x_mid(self):
        """
        Average x of rod.
        """
        return self._x_mid

    @property
    def y_mid(self):
        """
        Average y of rod.
        """
        return self._y_mid

    @property
    def kappa(self):
        """
        L/D of rod.
        """
        return self._kappa

    @property
    def angle(self):
        """
        Angle of rod.
        """
        return self._angle

    def is_in_circle(self, center, rad):
        """
        Checks if rod is in the circle defined by the given center and
        the given rad.
        """
        return is_in_circle(self.x_mid, self.y_mid,
                            center[0], center[1], rad)

    def has_valid_proportions(self, kappas, allowed_error):
        """
        Checks if rod has valid L/D (kappas are possibles values
        for L/D).
        """
        passed = []
        try:
            for kappa in kappas:
                condition = abs(self.kappa-kappa) < allowed_error
                passed.append(condition)
#kappa is not an array, so there is only 1 kappa.
        except TypeError:
            condition = abs(self.kappa-kappas) < allowed_error
            passed.append(condition)
        output = False
        for condition in passed:
            output = output or condition
        return output

    def is_valid_rod(self, kappas,
                    allowed_kappa_error,
                    zone_coords):
        """
        Check if rod is valid checking L/D and distance to center.
        TODO: If rods are near, kappa is not correct.
        """
        center = (zone_coords[0], zone_coords[1])
        is_in_main = self.is_in_circle(center, zone_coords[2])
        has_valid_proportions = self.has_valid_proportions(kappas,
                                                           allowed_kappa_error)
        output = is_in_main and has_valid_proportions
        return output

    def distance_to_rod(self, rod):
        """
        Returns the distance to another rod.
        """
        diff_x = abs(self.x_mid-rod.x_mid)
        diff_y = abs(self.y_mid-rod.y_mid)
        return math.sqrt(diff_x**2+diff_y**2)

    def angle_between_rods(self, rod):
        """
        Returns value of angle that formes this rod with another.
        """
        angle1 = (self.angle-rod.angle)*math.pi/180
        angle2 = math.pi - (self.angle-rod.angle)*math.pi/180
        return min(angle1, angle2)

