set pm3d
set hidden3d
set term png
set output "12r.png"
splot "12-epa.data"  u 1:2:4:4 w l t "real"

set output "12i.png"
splot "12-epa.data"  u 1:2:5:5 w l t "imaginary"

set output "12a.png"
splot "12-epa.data"  u 1:2:6 w l t "absolute"
