import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sqlite3 as sql
from matplotlib import animation

resol = 10
sigma = 1

connection = sql.connect("rods.db")
cursor = connection.cursor()
cursor2 = connection.cursor()
get_rods_sql = "select xmid,ymid,major,minor from datos where experiment_id=? and file_id=?"
get_file_ids_sql = "select file_id from datos where experiment_id=?"

cursor.execute("select distinct experiment_id from datos")
experiment_ids = cursor.fetchall()
experiment_ids = [experiment_id[0] for experiment_id in experiment_ids]


for experiment_id in experiment_ids:
    file_ids = cursor.execute(get_file_ids_sql, (experiment_id,))
    for file_id in file_ids:
        distr = np.array([[0.0 for dummy1 in range(resol)] for dummy2 in range(resol)])
        cursor2.execute(get_rods_sql, (experiment_id, file_id[0]))
        rods = cursor2.fetchall()
        rods = np.array([np.array(rod) for rod in rods])
        num_rods = len(rods)
        xvals = rods[:,0]
        yvals = rods[:,1]
        x0 = min(xvals)
        y0 = min(yvals)
        dx = float(max(xvals)-x0)/resol
        dy = float(max(yvals)-y0)/resol
        rad = (dx*resol/2.0 + dy*resol/2.0)/2.0
        center = np.array([(min(xvals)+max(xvals))/2.0, (min(yvals)+max(yvals))/2.0])
        for row in range(resol):
            for col in range(resol):
                pos = np.array([col*dx+x0, row*dy+y0])
                if np.sum((pos-center)**2) <= rad**2:
                    val = -((xvals-col*dx-x0)**2/dx**2+np.abs(yvals-row*dy-y0)**2/dy**2)
                    exps = np.exp(val/(4*sigma**2))/np.sqrt(2*math.pi)*sigma
                    distr[row,col] += np.sum(exps)/num_rods
        x_vals = np.array(range(resol))
        y_vals = np.array(range(resol))
        mesh = np.array(np.meshgrid(x_vals, y_vals))
        x_vals, y_vals = tuple(mesh.reshape(2,resol**2))
        distr = list(distr.reshape(1,len(distr)**2)[0])
        size = 400
        plt.scatter(x_vals, y_vals, c=distr, marker='s', s=size)
        plt.show()

connection.commit()
connection.close()
