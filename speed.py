"""
Speed animation creation module.
"""

import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sqlite3 as sql
from matplotlib import animation
import matplotlib
from multiprocessing import Pool
import itertools

np.warnings.filterwarnings('ignore')

_writer = animation.writers['ffmpeg']
writer = _writer(fps=15, metadata=dict(artist='Jaime Perez Aparicio'), bitrate=1800)
WRITER = writer
CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
REWRITE_LAST = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

resol = 30
cresol = 256*32
sigma = 5

connection = sql.connect("rods.db")
cursor = connection.cursor()
cursor2 = connection.cursor()

get_rods_sql = "select xmid,ymid,major/minor,angle from datos where experiment_id=? and file_id=? order by ymid,xmid"
get_file_ids_sql = "select distinct file_id from datos where experiment_id=?"
get_close_rods = "select xmid, ymid, angle from datos where experiment_id=? "
get_close_rods+= "and file_id=? and ymid between ? and ? and xmid between ? and ? and angle between ? and ? "
get_close_rods+= "and major/minor between ? and ? order by ((xmid-?)*(xmid-?)+(ymid-?)*(ymid-?))*angle" 

cursor.execute("select distinct experiment_id from datos")
experiment_ids = cursor.fetchall()
experiment_ids = [experiment_id[0] for experiment_id in experiment_ids]

colors = plt.cm.jet(np.linspace(-1,1,cresol))
white = (1,1,1,1)
plot_ = None
fig = plt.figure()

dd = 30
dangle = 10
dk = 6

def _speeds_matrix(experiment_id_, file_ids, index, x_vals, y_vals):
    """
    Computes linear_speed+2*angular_speed/3 for each grid square.
    """
    cursor2.execute(get_rods_sql, (str(experiment_id_), str(file_ids[index])))
    rods0 = np.array(cursor2.fetchall())
    xmin = min(rods0[:, 0])
    xmax = max(rods0[:, 0])
    ymin = min(rods0[:, 1])
    ymax = max(rods0[:, 1])
    rad = (float(xmax-xmin)+float(ymax-ymin))/4.
    center = np.array([float(xmax+xmin)/2., float(ymax+ymin)/2.])
    dx = float(xmax-xmin)/resol
    dy = float(ymax-ymin)/resol
    speeds = [[np.nan for dummy in range(resol)] for dummy2 in range(resol)]
    empties = 0.0
    for rod in rods0:
        pos0 = np.array([rod[0], rod[1]])
        angle = rod[3]
        k0 = rod[2]
        args = (experiment_id_, file_ids[index+1], pos0[1]-dd,
                pos0[1]+dd, pos0[0]-dd, pos0[0]+dd,
                angle-dangle, angle+dangle, k0-dk, k0+dk,
                pos0[0], pos0[0], pos0[1], pos0[1])
        cursor2.execute(get_close_rods, args)
        rodsf = np.array(cursor2.fetchall())
        try:
            rodf = rodsf[0]
            dist2 = np.sum((rodf[:2]-rod[:2])**2)
            speed = np.sqrt(dist2)
            ang_speed = np.radians(abs(rodf[2]-rod[2]))
            ang_speed = min(ang_speed, np.pi-ang_speed)
            ang_speed = min(ang_speed, np.pi/2-ang_speed)
            for index_x in range(resol):
                for index_y in range(resol):
                    xx = index_x*dx+xmin
                    yy = index_y*dx+ymin
                    rr = np.array([xx, yy])
                    dist = np.sqrt(np.sum((rr-center)**2))
                    if dist<=rad:
                        ddx = rod[0]-xx
                        ddy = rod[1]-yy
                        ddr = np.array([ddx, ddy])
                        value = np.exp(-(np.sum(ddr**2))/(2.*dx**2*sigma**2))/(np.sqrt(2*np.pi)*sigma)*(speed+2.0*ang_speed/3)
                        if np.isnan(speeds[index_x][index_y]):
                            speeds[index_x][index_y] = value
                        else:
                            speeds[index_x][index_y] += value
        except IndexError:
            pass
    return np.array(speeds)
    

def _get_distr(experiment_id_, file_ids, x_vals, y_vals):
    """
    Computes plottable data average.
    """
    cursor2 = connection.cursor()
    speeds = []
    for index in range(4):
        file_id = file_ids[index]
        speed = _speeds_matrix(experiment_id_, file_ids, index, x_vals, y_vals)
        speeds.append(speed)
    return np.array(sum(speeds)/4.0)

def _plot_frame(experiment_id_, file_ids, plot, fig):
    """
    Plots frame.
    """
    plt.clf()
    x_vals, y_vals = np.array(range(resol)), np.array(range(resol))
    distr = _get_distr(experiment_id_, file_ids, x_vals, y_vals)
    #first plot. It create axes...
    mesh = np.array(np.meshgrid(x_vals, y_vals))
    x_vals, y_vals = tuple(mesh.reshape(2, resol**2))
    size = 2000.0/resol
    elements = distr.shape[0]*distr.shape[1]
    distr.reshape((1, elements))
    plot = plt.scatter(x_vals, y_vals, c=distr, marker='s', s=size,
                       vmin=0, vmax=2)
    cb = plt.colorbar()
    cb.set_label('density change')
    plt.xlabel("x [norm]")
    plt.ylabel("y [norm]")
    plt.suptitle('speeds distribution')
    return plot

def create_animation(experiment_id):
    """
    Creates animation for an experiment.
    """
    file_ids = cursor.execute(get_file_ids_sql, (experiment_id,)).fetchall()
    file_ids = [file_id[0] for file_id in file_ids]
    num_frames = int(len(file_ids)/5)
    global plot_
    global fig
    plot_ = None
    exit = False
    frame_idx = 0
    msg = "Computing speeds matrix for experiment {}".format(experiment_id)
    print(msg)
    def animate(frame_idx):
        """
        Wrapper
        """
        global plot_
        global fig
        progress = str(frame_idx) + "/" + str(num_frames) + "\t("
        progress += "{0:.2f}%)"
        print(progress.format(frame_idx*100.0/num_frames),end='\r')
        file_id0 = frame_idx*5
        plot_ = _plot_frame(int(experiment_id), file_ids[file_id0:file_id0+5] , plot_, fig)
    anim = animation.FuncAnimation(fig, animate, frames=num_frames, repeat=False)
    name = 'speeds_animation{}.mp4'.format(experiment_id)
    try:
        anim.save(name)
    except KeyboardInterrupt:
        connection.commit()
        connection.close()
        raise KeyboardInterrupt
    connection.commit()

for experiment_id in experiment_ids:
    create_animation(experiment_id)

connection.close()
