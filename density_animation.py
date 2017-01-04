import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sqlite3 as sql
from matplotlib import animation
import matplotlib


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


def _plot_frame(experiment_id_, file_ids, plot):
    """
    Plots a frame doing some averages over burst.
    """
    distr = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
    for index in range(5):
        distr_6 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
        distr_12 = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
        cursor2.execute(get_rods_sql, (str(experiment_id_), str(file_ids[index])))
        rods = cursor2.fetchall()
        #len(file_ids) can be < 5
        if not len(rods):
            print("Empty")
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
        for row in range(resol):
            for col in range(resol):
                pos = np.array([col*dx+x0, row*dy+y0])
                if np.sum((pos-center)**2) <= rad**2:
                    val_6 = -((rods_6[:, 0]-col*dx-x0)**2/dx**2+np.abs(rods_6[:, 1]-row*dy-y0)**2/dy**2)
                    exps_6 = np.exp(val_6/(4*sigma**2))/np.sqrt(2*math.pi)*sigma
                    distr_6[row, col] += np.sum(exps_6)/num_rods
                    val_12 = -((rods_12[:, 0]-col*dx-x0)**2/dx**2+np.abs(rods_12[:, 1]-row*dy-y0)**2/dy**2)
                    exps_12 = np.exp(val_12/(4*sigma**2))/np.sqrt(2*math.pi)*sigma
                    distr_12[row, col] += np.sum(exps_12)/num_rods
        #average over burst
        distr += (distr_12-distr_6)/(distr_12+distr_6)/5.0
    distr = list(distr.reshape(1, len(distr)**2)[0])
    if not(plot is None):
        #only update colors.
        colors_ = [colors[int(value*cresol/2)] for value in distr]
        plot.set_color(colors)
        return plot
    else:
        #first plot. It create axes...
        x_vals, y_vals = np.array(range(resol)), np.array(range(resol))
        mesh = np.array(np.meshgrid(x_vals, y_vals))
        x_vals, y_vals = tuple(mesh.reshape(2, resol**2))
        size = 4000.0/resol
        plot = plt.scatter(x_vals, y_vals, c=distr, marker='s', s=size)
        cb = plt.colorbar()
        cb.set_label('(n_12-n_6)/(n_12+n_6)')
        plt.xlabel("x [norm]")
        plt.ylabel("y [norm]")
        plt.suptitle('density distribution')
        first = False
        return plot

def create_animation(experiment_id):
    """
    Creates animation for an experiment.
    """
    print(experiment_id)
    file_ids = cursor.execute(get_file_ids_sql, (experiment_id,)).fetchall()
    fig = plt.figure()
    num_frames = int(len(file_ids)/5)
    print("\n")
    plot = None
    def animate(frame_idx):
        """
        Wrapper
        """
        global plot
        progress = CLEAR_LAST + str(frame_idx*5) + "/" + str(num_frames) + "("
        progress += "{0:.2f}%)"
        print(progress.format(frame_idx*100.0/num_frames))
        file_id0 = frame_idx*5
        plot = _plot_frame(int(experiment_id), file_ids[file_id0:file_id0+5] , plot)
    anim = animation.FuncAnimation(fig, animate, frames=num_frames)
    name = 'density_animation{}.mp4'.format(experiment_id)
    anim.save(name)


for experiment_id in experiment_ids:
    create_animation(experiment_id)

connection.commit()
connection.close()
