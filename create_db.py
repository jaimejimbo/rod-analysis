import sqlite3 as sql

conn = sql.connection('rods.db')
cursor = conn.cursor()


cursor.execute("CREATE TABLE datos(experiment_id  integer, file_id integer, id integer, area real, xmid real, ymid real, major real, minor real, angle real,    feret real, feretx real, ferety real, feretangle real, minferet real, xstart real, ystart real);")
cursor.execute("CREATE TABLE times(experiment_id integer, file_id integer, time real);")
cursor.execute("CREATE INDEX index3 on times(experiment_id,file_id);")
cursor.execute("CREATE INDEX index1 on datos(experiment_id, file_id);")

conn.commit()
conn.close()
