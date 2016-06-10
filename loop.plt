
filename = "rods_".i
plotfile = "graph".i.".jpg"

print filename." ".plotfile
set output plotfile
set yrange [0:200]
plot filename using 1:8
set output
i=i+1
if (i <= n) reread
