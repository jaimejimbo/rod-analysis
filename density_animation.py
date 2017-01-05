"""
Density animation creation module.
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
CLEAR_LAST = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

resol = 60
cresol = 1000
sigma = 1

connection = sql.connect("rods.db")
cursor = connection.cursor()
cursor2 = connection.cursor()
get_rods_sql = "select xmid,ymid,major,minor from datos where experiment_id=? and file_id=?"
get_file_ids_sql = "select distinct file_id from datos where experiment_id=?"

cursor.execute("select distinct experiment_id from datos")
experiment_ids = cursor.fetchall()
experiment_ids = [experiment_id[0] for experiment_id in experiment_ids]

colors = plt.cm.jet(np.linspace(-1,1,cresol))
plot_ = None
fig = plt.figure()


def pool_func(row, col, rods_6, rods_12, dx, dy, x0, y0, center, rad, num_rods):
    """
    Function to be used in pool.
    """
    pos = np.array([col*dx+x0, row*dy+y0])
    if np.sum((pos-center)**2) <= rad**2:
        distr_6 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
        distr_12 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
        val_6 = -((rods_6[:, 0]-col*dx-x0)**2/dx**2+np.abs(rods_6[:, 1]-row*dy-y0)**2/dy**2)
        exps_6 = np.exp(val_6/(4*sigma**2))/np.sqrt(2*math.pi)*sigma
        distr_6[row, col] += np.sum(exps_6)/num_rods
        val_12 = -((rods_12[:, 0]-col*dx-x0)**2/dx**2+np.abs(rods_12[:, 1]-row*dy-y0)**2/dy**2)
        exps_12 = np.exp(val_12/(4*sigma**2))/np.sqrt(2*math.pi)*sigma
        distr_12[row, col] += np.sum(exps_12)/num_rods
        #average over burst
        return (distr_12-distr_6)/(distr_12+distr_6)/5.0
    return np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    
def _get_distr(experiment_id_, file_ids):
    """
    Computes plottable data.
    """
    distr = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    for index in range(5):
        cursor2.execute(get_rods_sql, (str(experiment_id_), str(file_ids[index])))
        rods = cursor2.fetchall()
        #len(file_ids) can be < 5
        if not len(rods):
            print("Empty")
            print(experiment_id_)
            print(file_ids)
            return
        rods = np.array([np.array(rod) for rod in rods])
        rods_6, rods_12 = [], []
        for rod in rods:
            if 6 < float(rod[2])/rod[3] < 12:
                rods_6.append(rod)
            elif 12 < float(rod[2])/rod[3] < 20:
                rods_12.append(rod)
        rods_6, rods_12 = np.array(rods_6), np.array(rods_12)
        num_rods = len(rods)
        xvals, yvals = rods[:, 0], rods[:, 1]
        x0, y0 = min(xvals), min(yvals)
        dx = float(max(xvals)-x0)/resol
        dy = float(max(yvals)-y0)/resol
        rad = (dx*resol/2.0 + dy*resol/2.0)/2.0
        center = np.array([(min(xvals)+max(xvals))/2.0, (min(yvals)+max(yvals))/2.0])
        rows, cols = np.meshgrid(range(resol), range(resol))
        rows = np.array(rows).reshape(1,resol**2)[0]
        cols = np.array(cols).reshape(1,resol**2)[0]
        pool = Pool(4)
        r = itertools.repeat
        izip = zip(rows, cols, r(rods_6),
                   r(rods_12), r(dx), r(dy),
                   r(x0), r(y0), r(center),
                   r(rad), r(num_rods))
        distr = sum(pool.starmap(pool_func, izip))
    print(distr)
#all values are nan
    distr[np.argwhere(np.isnan(distr))] = -2
    distr = list(distr.reshape(1, len(distr)**2)[0])
    return distr
    

def _plot_initial_frame(experiment_id_, file_ids, plot, fig):
    """
    Plots initial frame.
    """
    distr = _get_distr(experiment_id_, file_ids)
    #first plot. It create axes...
    x_vals, y_vals = np.array(range(resol)), np.array(range(resol))
    mesh = np.array(np.meshgrid(x_vals, y_vals))
    x_vals, y_vals = tuple(mesh.reshape(2, resol**2))
    size = 2000.0/resol
    plot = plt.scatter(x_vals, y_vals, c=distr, marker='s', s=size)
    #cb = plt.colorbar()
    #cb.set_label('(n_12-n_6)/(n_12+n_6)')
    plt.xlabel("x [norm]")
    plt.ylabel("y [norm]")
    plt.suptitle('density distribution')
    #fig.draw()
    #plt.show()
    return plot

def _plot_frame(experiment_id_, file_ids, plot, fig):
    """
    Updates the frame.
    """
    #only update colors.
    distr = _get_distr(experiment_id_, file_ids)
    colors_ = [colors[int(value*cresol/2)] for value in distr]
    plot.set_color(colors)
    #fig.draw()
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
    file_id0 = frame_idx*5
    plot_ = _plot_initial_frame(int(experiment_id), file_ids[file_id0:file_id0+5] , plot_, fig)
    def animate(frame_idx):
        """
        Wrapper
        """
        global plot_
        global fig
        progress = CLEAR_LAST + str(frame_idx) + "/" + str(num_frames) + "("
        progress += "{0:.2f}%)"
        print(progress.format(frame_idx*100.0/num_frames))
        file_id0 = frame_idx*5
        if plot_ is None:
            plot_ = _plot_initial_frame(int(experiment_id), file_ids[file_id0:file_id0+5] , plot_, fig)
        else:
            plot_ = _plot_frame(int(experiment_id), file_ids[file_id0:file_id0+5] , plot_, fig)
    try:
        anim = animation.FuncAnimation(fig, animate, frames=num_frames, repeat=False)
    except KeyboardInterrupt:
        exit = True
    name = 'density_animation{}.mp4'.format(experiment_id)
    anim.save(name)
    if exit:
        raise KeyboardInterrupt

for experiment_id in experiment_ids:
    create_animation(experiment_id)

connection.commit()
connection.close()
