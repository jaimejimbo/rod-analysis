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
REWRITE_LAST = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

resol = 50
cresol = 256*32
sigma = 0.5

connection = sql.connect("rods.db")
cursor = connection.cursor()
cursor2 = connection.cursor()
get_rods_sql = "select xmid,ymid,major,minor from datos where experiment_id=? and file_id=? order by ymid,xmid"
get_file_ids_sql = "select distinct file_id from datos where experiment_id=?"

cursor.execute("select distinct experiment_id from datos")
experiment_ids = cursor.fetchall()
experiment_ids = [experiment_id[0] for experiment_id in experiment_ids]

colors = plt.cm.jet(np.linspace(-1,1,cresol))
white = (1,1,1,1)
plot_ = None
fig = plt.figure()


def _update_matrix(row, col, rods_6, rods_12, dx, dy, x0, y0, center, rad, num_rods, distr6, distr12):
    """
    Updates matrix with density values.
    """
    pos = np.array([col*dx+x0, row*dy+y0])
    if np.sum((pos-center)**2) <= rad**2:
        xcenter = col*dx-x0
        ycenter = row*dy-y0
        xvals6 = rods_6[:, 0]
        yvals6 = rods_6[:, 1]
        xvals12 = rods_12[:, 0]
        yvals12 = rods_12[:, 1]
        val6 = -((xvals6-xcenter)**2/dx**2 + (-ycenter)**2/dy**2)
        exps6 = np.exp(val6/(4*sigma**2))/(np.sqrt(2*math.pi)*sigma)
        distr6[row, col] = np.sum(exps6)/num_rods
        val12 = -((xvals12-xcenter)**2/dx**2 + (yvals12-ycenter)**2/dy**2)
        exps12 = np.exp(val12/(4*sigma**2))/(np.sqrt(2*math.pi)*sigma)
        distr12[row, col] = np.sum(exps12)/num_rods
    return distr6, distr12

def _density_matrix_process(experiment_id_, file_ids, rods):
    """
    Density matrix computing process.
    """
    distr_6 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    distr_12 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
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
            rods_6.append(rod[:2])
        elif 12 < float(rod[2])/rod[3] < 20:
            rods_12.append(rod[:2])
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
    for row, col in zip(rows, cols):
        distr_6, distr_12 = _update_matrix(row, col, rods_6, rods_12,
                      dx, dy, x0, y0, center, rad, num_rods,
                      distr_6, distr_12)
    return distr_6/5.0, distr_12/5.0

def _get_distr(experiment_id_, file_ids):
    """
    Computes plottable data.
    """
    cursor2 = connection.cursor()
    rods = []
    for index in range(5):
        cursor2.execute(get_rods_sql, (str(experiment_id_), str(file_ids[index])))
        rods.append(cursor2.fetchall())
    try:
        pool = Pool(4)
        r = itertools.repeat
        izip = zip(r(experiment_id_), r(file_ids), rods)
        distrs = np.array(pool.starmap(_density_matrix_process, izip))
    finally:
        pool.close()
        pool.join()
    distr_6 = sum(distrs[:,0,:,:])
    distr_12 = sum(distrs[:,1,:,:])
    distr = (distr_12-distr_6)/(distr_12+distr_6)
    distr = list(distr.reshape(1, len(distr)**2)[0])
    return distr

def _create_colors_proc(value):
    """
    Process
    """
    if np.isnan(value):
        color = white
    else:
        color = colors[int((value+1)*cresol/2)-1]
    return color

def _create_colors(distr):
    """
    Creates a color distribution for data.
    """
    colors_ = []
    try:
        pool = Pool(4)
        colors_ = pool.map(_create_colors_proc, distr)
    finally:
        pool.close()
        pool.join()
    return colors_

def _plot_frame(experiment_id_, file_ids, plot, fig):
    """
    Plots frame.
    """
    plt.clf()
    distr = _get_distr(experiment_id_, file_ids)
    #first plot. It create axes...
    x_vals, y_vals = np.array(range(resol)), np.array(range(resol))
    mesh = np.array(np.meshgrid(x_vals, y_vals))
    x_vals, y_vals = tuple(mesh.reshape(2, resol**2))
    #x_vals, y_vals, distr = delete_nones(x_vals, y_vals, distr)
    size = 2000.0/resol
    plot = plt.scatter(x_vals, y_vals, c=distr, marker='s', s=size,
                       vmin=-1, vmax=1)
    colors_ = _create_colors(distr)
    plot.set_color(colors_)
    cb = plt.colorbar()
    cb.set_label('(n_12-n_6)/(n_12+n_6)')
    plt.xlabel("x [norm]")
    plt.ylabel("y [norm]")
    plt.suptitle('density distribution')
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
    msg = "Computing densiry matrix for experiment {}".format(experiment_id)
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
    name = 'density_animation{}.mp4'.format(experiment_id)
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
