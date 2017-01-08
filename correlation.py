"""
Correlation animation creation module.
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
sigma = 1.0

connection = sql.connect("rods.db")
cursor = connection.cursor()
cursor2 = connection.cursor()
get_rods_sql = "select xmid,ymid,major,minor,angle from datos where experiment_id=? and file_id=? order by ymid,xmid"
get_file_ids_sql = "select distinct file_id from datos where experiment_id=?"

cursor.execute("select distinct experiment_id from datos")
experiment_ids = cursor.fetchall()
experiment_ids = [experiment_id[0] for experiment_id in experiment_ids]

colors = plt.cm.jet(np.linspace(-1,1,cresol))
white = (1,1,1,1)
fig, axarr = plt.subplots(2, 2, sharex=True, sharey=True)


def _update_matrix(row, col, rods_6, rods_12, dx, dy, x0, y0, center, rad, num_rods, distr6_q2, distr6_q4, distr12_q2, distr12_q4):
    """
    Updates matrix with density values.
    """
    pos = np.array([col*dx+x0, row*dy+y0])
    if np.sum((pos-center)**2) <= rad**2:
        xcenter = col*dx+x0
        ycenter = row*dy+y0
        xvals6 = rods_6[:, 0]
        yvals6 = rods_6[:, 1]
        xvals12 = rods_12[:, 0]
        yvals12 = rods_12[:, 1]
        angles6 = rods_6[:, 2]
        angles12 = rods_12[:, 2]
        q2_6 = (np.cos(2*angles6) + np.sin(2*angles6))/num_rods
        q2_12 = (np.cos(2*angles12) + np.sin(2*angles12))/num_rods
        q4_6 = (np.cos(4*angles6) + np.sin(2*angles6))/num_rods
        q4_12 = (np.cos(4*angles12) + np.sin(4*angles12))/num_rods
        val6 = -((xvals6-xcenter)**2/dx**2 + (yvals6-ycenter)**2/dy**2)
        exps6 = np.exp(val6/(4*sigma**2))/(np.sqrt(2*math.pi)*sigma)
        distr6_q2[row, col] = np.dot(q2_6, exps6)/5.0
        distr6_q4[row, col] = np.dot(q4_6, exps6)/5.0
        val12 = -((xvals12-xcenter)**2/dx**2 + (yvals12-ycenter)**2/dy**2)
        exps12 = np.exp(val12/(4*sigma**2))/(np.sqrt(2*math.pi)*sigma)
        distr12_q2[row, col] = np.dot(q2_12, exps12)/5.0
        distr12_q4[row, col] = np.dot(q4_12, exps12)/5.0
    else:
        distr6_q2[row, col] = np.nan
        distr6_q4[row, col] = np.nan
        distr12_q2[row, col] = np.nan
        distr12_q4[row, col] = np.nan
    return distr6_q2, distr6_q4, distr12_q2, distr12_q4

def _correlation_matrix_process(experiment_id_, file_ids, rods):
    """
    Order parameter matrix computing process.
    """
    distr6_q2 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    distr6_q4 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    distr12_q2 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    distr12_q4 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
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
            rods_6.append(np.append(rod[:2],rods[4]))
        elif 12 < float(rod[2])/rod[3] < 20:
            rods_12.append(np.append(rod[:2],rods[4]))
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
        distr6_q2, distr6_q4, distr12_q2, distr12_q4 = _update_matrix(row, col, rods_6, 
                      rods_12, dx, dy, x0, y0, center, rad, num_rods,
                      distr6_q2, distr6_q4, distr12_q2, distr12_q4)
    return distr6_q2, distr6_q4, distr12_q2, distr12_q4

def _get_distrs(experiment_id_, file_ids):
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
        distrs = np.array(pool.starmap(_correlation_matrix_process, izip))
    finally:
        pool.close()
        pool.join()
    distr6_q2 = sum(distrs[:,0,:,:])
    distr6_q4 = sum(distrs[:,1,:,:])
    distr12_q2 = sum(distrs[:,2,:,:])
    distr12_q4 = sum(distrs[:,3,:,:])
    return distr6_q2, distr6_q4, distr12_q2, distr12_q4

def _plot_frame(experiment_id_, file_ids, fig, axarr):
    """
    Plots frame.
    """
    distr6_q2, distr6_q4, distr12_q2, distr12_q4 = _get_distrs(experiment_id_, file_ids)
    x_vals, y_vals = np.array(range(resol)), np.array(range(resol))
    mesh = np.array(np.meshgrid(x_vals, y_vals))
    x_vals, y_vals = tuple(mesh.reshape(2, resol**2))
    rad = (max(x_vals)-min(x_vals))/(2.0*resol)
    size = 125*(rad)**2
    sp1 = axarr[0, 0].scatter(x_vals, y_vals, c=distr6_q2, marker='s', s=size)
    sp2 = axarr[0, 1].scatter(x_vals, y_vals, c=distr6_q4, marker='s', s=size)
    sp3 = axarr[1, 0].scatter(x_vals, y_vals, c=distr12_q2, marker='s', s=size)
    sp4 = axarr[1, 1].scatter(x_vals, y_vals, c=distr12_q4, marker='s', s=size)
    fig.suptitle("angular correlations")
    axarr[0, 0].set_title("K6 Q2")
    axarr[0, 1].set_title("K6 Q4")
    axarr[1, 0].set_title("K12 Q2")
    axarr[1, 1].set_title("K12 Q4")
    for ax_ in axarr.reshape(1, 4):
        for ax in ax_:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
    cb1 = plt.colorbar(sp1, cax=axarr[0, 0])
    cb2 = plt.colorbar(sp2, cax=axarr[0, 1])
    cb3 = plt.colorbar(sp3, cax=axarr[1, 0])
    cb4 = plt.colorbar(sp4, cax=axarr[1, 1])
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)    

def create_animation(experiment_id):
    """
    Creates animation for an experiment.
    """
    file_ids = cursor.execute(get_file_ids_sql, (experiment_id,)).fetchall()
    file_ids = [file_id[0] for file_id in file_ids]
    num_frames = int(len(file_ids)/5)
    global fig
    global axarr
    exit = False
    frame_idx = 0
    msg = "Computing order parameter matrix for experiment {}".format(experiment_id)
    print(msg)
    def animate(frame_idx):
        """
        Wrapper
        """
        global fig
        global axarr
        progress = str(frame_idx) + "/" + str(num_frames) + "\t("
        progress += "{0:.2f}%)"
        print(progress.format(frame_idx*100.0/num_frames),end='\r')
        file_id0 = frame_idx*5
        _plot_frame(int(experiment_id), file_ids[file_id0:file_id0+5], fig, axarr)
    anim = animation.FuncAnimation(fig, animate, frames=num_frames, repeat=False)
    name = 'correlation_animation{}.mp4'.format(experiment_id)
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
