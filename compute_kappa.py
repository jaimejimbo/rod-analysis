import sqlite3 as sql
import numpy as np

conn = sql.connect("./rods.db")
cursor1 = conn.cursor()
cursor1.execute("select experiment_id,file_id,id,major,minor from datos")
cursor2 = conn.cursor()

while True:
    data = cursor1.fetchone()
    if data is None:
        break
    [experiment_id, file_id, rod_id, major, minor] = data
    kappa = float(major)/minor
    data = (experiment_id,file_id,rod_id,kappa)
    cursor2.execute("insert into kappas values (?,?,?,?)",data)

conn.commit()
conn.close()
    
