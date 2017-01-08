"""
Rods' animation creation module.
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
get_rods_sql = "select xmid,ymid,major,minor,angle from datos where experiment_id=? and file_id=? order by ymid,xmid"
get_file_ids_sql = "select distinct file_id from datos where experiment_id=?"

cursor.execute("select distinct experiment_id from datos")
experiment_ids = cursor.fetchall()
experiment_ids = [experiment_id[0] for experiment_id in experiment_ids]

#colors = plt.cm.jet(np.linspace(-1,1,cresol))
#white = (1,1,1,1)
plot12, plot6 = None, None
fig = plt.figure()
initial = True

def _plot_frame(experiment_id_, file_id, plot12, plot6, fig):
    """
    Plots frame.
    """
    plt.clf()
    rods = cursor.execute(get_rods_sql, (str(experiment_id_), str(file_id))).fetchall()
    rods = np.array(rods)
    conds = rods[:, 2]/rods[:, 3]>=12
    lengths = np.where(conds, 12, 6)
    scale = np.average(rods[:, 2]/lengths)/2
    lengths = lengths*scale
    dx = lengths*np.cos(-np.radians(rods[:, 4]))
    dy = lengths*np.sin(-np.radians(rods[:, 4]))
    pd = np.array([[rods[:, 0]-dx, rods[:, 1]-dy], [rods[:, 0]+dx, rods[:, 1]+dy]])
    pd12 = pd[:, :, np.where(conds)[0]]
    pd6 = pd[:, :, np.where(np.logical_not(conds))[0]]
    plot12 = plt.plot([pd12[0][0], pd12[1][0]], [pd12[0][1], pd12[1][1]], color='r', linewidth=2.5)
    plot6 = plt.plot([pd6[0][0], pd6[1][0]], [pd6[0][1], pd6[1][1]], color='b', linewidth=2.5)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.suptitle('rods distribution')
    return plot12, plot6

def create_animation(experiment_id):
    """
    Creates animation for an experiment.
    """
    file_ids = cursor.execute(get_file_ids_sql, (experiment_id,)).fetchall()
    file_ids = [file_id[0] for file_id in file_ids]
    num_frames = len(file_ids)
    global plot12, plot6
    global fig
    global initial
    plot_ = None
    exit = False
    initial = True
    frame_idx = 0
    msg = "Plotting rods for experiment {}".format(experiment_id)
    print(msg)
    def animate(frame_idx):
        """
        Wrapper
        """
        global plot12, plot6
        global fig
        global initial
        progress = str(frame_idx) + "/" + str(num_frames) + "\t("
        progress += "{0:.2f}%)"
        print(progress.format(frame_idx*100.0/num_frames),end='\r')
        plot12, plot6 = _plot_frame(int(experiment_id), file_ids[frame_idx] , plot12, plot6, fig)
    anim = animation.FuncAnimation(fig, animate, frames=num_frames, repeat=False)
    name = 'rods_animation{}.mp4'.format(experiment_id)
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
