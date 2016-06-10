set term png enhanced size 1200,1200
set output 'seg.png'
set size ratio 1
p [400:2200][1800:0] 'clusters_long_5640'  u 1:2 w l lw 4 lc 23 t '', 'clusters_short_5640'  u 1:2 w l lw 4 lc 7 t ''
