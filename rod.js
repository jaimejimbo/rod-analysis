/*	Rod object. See code for more info.
 *  Author: Jaime Pérez Aparicio
 *  email: 19.jaime.91@gmail.com
 *  license: GPL
 */ 

var rod = function(ID, area, xm, ym, major, minor, angle, feret, feretx, ferety, feretangle, minferet, xstart, ystart){
    this.ID = ID;
    this.area = area;
    this.xm = xm;
    this.ym = ym;
    this.major = major;
    this.minor = minor;
    this.angle = angle;
    this.feret = feret;
    this.feretx = feretx;
    this.ferety = ferety;
    this.feretangle = feretangle;
    this.minferet = minferet;
    this.xstart = xstart;
    this.ystart = ystart;
    this.LD = this.major/this.minor;
};

