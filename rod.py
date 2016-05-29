"""
Rod object.
"""
import math, methods, matrix

class Rod(object):
    """
    Rod object.
    """

    def __init__(self, args_tuple, kappa=None, real_length=None):
        """
        Initialization of rod
        """
        (ID, area, xmid, ymid, major, minor,
                        angle, feret, feretx, ferety,
                        feretangle, minferet, xstart, ystart) = args_tuple
                                                #Column
        self._id = int(ID)                      #0
        self._area = float(area)                #1
        self._x_mid = float(xmid)               #2
        self._y_mid = float(ymid)               #3
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
        self._direction_matrix = matrix.zeros(2, 2)
        self._kappa = float(feret)/float(minferet)
        self._real_kappa = kappa
        self._real_length = real_length
        #self._feret = float(real_length)
        #self._min_feret = real_length/self._kappa

    @property
    def area(self):
        """
        Returns area covered by rod.
        """
        return float(self._real_length**2)/self._real_kappa

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

    @property
    def center(self):
        """
        Returns position of the center of the rod.
        """
        return self._x_mid, self._y_mid

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

    def __repr__(self):
        """
        String transformation.
        """
        output = ""
        output += "id: "+str(self.identifier)+"\n"
        output += "center: "+str(self.center)+"\n"
        output += "angle: "+str(self.angle)+"\n"
        return output

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
        Computed L/D of rod.
        """
        return self._kappa

    @property
    def kappa(self):
        """
        Real L/D of rod. 
        """
        return self._kappa

    @kappa.setter
    def kappa(self, value):
        """
        L/D of rod.
        """
        self._kappa = value

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
        return methods.is_in_circle(self.x_mid, self.y_mid,
                            center[0], center[1], rad)

    def has_valid_proportions(self, kappas, allowed_error):
        """
        Checks if rod has valid L/D (kappas are possibles values
        for L/D).
        """
        try:
            for kappa in kappas:
                if abs(self.kappa-kappa) < allowed_error:
                    return True
            return False
        except TypeError:
            return abs(self.kappa-kappas) < allowed_error

    def is_valid_rod(self, kappas,
                    allowed_kappa_error,
                    zone_coords):
        """
        Check if rod is valid checking L/D and distance to center.
        """
        center = (zone_coords[0], zone_coords[1])
        is_in_main = self.is_in_circle(center, zone_coords[2])
        has_valid_proportions = self.has_valid_proportions(kappas,
                                                           allowed_kappa_error)
        output = is_in_main and has_valid_proportions
        if output:
            self._kappa = self._real_kappa
            self._feret = self._real_length
            self._min_feret = self._feret*1.0/self._kappa
        return output

    def is_valid_rod_length(self, length,
                        length_error, kappa,
                        zone_coords):
        """
            Check if rod's length is correct and if it's in the circle.
        """
        center = (zone_coords[0], zone_coords[1])
        is_in_main = self.is_in_circle(center, zone_coords[2])
        rod_length = self.feret
        valid_length = ((length-length_error) <= rod_length <= (length+length_error))
        cond = valid_length and is_in_main
        if cond:
            self._kappa = kappa
            self._feret_min = self._feret/float(kappa)
            self._feret = self._real_length
        return cond

    def vector_to_rod(self, rod, scale=1):
        """
        Returns a vector that joins 2 rods.
        Start in self.
        """
        scale = float(scale)
        diff_x = (rod.x_mid-self.x_mid)*scale
        diff_y = (rod.y_mid-self.y_mid)*scale
        return (diff_x, diff_y)

    def distance_to_rod(self, rod, scale=1):
        """
        Returns the distance to another rod.
        """
        return methods.vector_module(self.vector_to_rod(rod, scale=scale))

    def angle_between_rods(self, rod):
        """
        Returns value of angle that formes this rod with another.
        """
        angle = abs(self.angle-rod.angle) % 180
        return angle

    @property
    def direction_matrix(self):
        """
        Returns a matrix with the form:
        ex^2-1  ex*ey
        ex*ey   ey^2-1
        """
        if self._direction_matrix == matrix.zeros(2, 2):
            e_x = math.cos(self.angle)
            e_y = math.sin(self.angle)
            self._direction_matrix[0][0] = 2*e_x**2-1
            self._direction_matrix[1][1] = 2*e_y**2-1
            self._direction_matrix[1][0] = e_x*e_y
            self._direction_matrix[0][1] = e_x*e_y
        return self._direction_matrix

    @property
    def real_length(self):
        """
        Length of rod (pixels)
        """
        return self._real_length

