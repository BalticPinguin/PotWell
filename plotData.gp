#set pm3d
#splot "out2.data"  u 1:2:4:4 title "real"
#pause -1
#
#splot "out2.data"  u 1:2:5:5 title "imaginary"
#pause -1
#
#splot "out2.data"  u 1:2:6 title "absolute"
#pause -1

# set pm3d
# splot "out2a.data"  u 1:2:4:4 title "real"
# pause -1
# 
# splot "out2a.data"  u 1:2:5:5 title "imaginary"
# pause -1
# 
# splot "out2a.data"  u 1:2:6 title "absolute"
# pause -1
#set pm3d
#splot "0_ep.data"  u 1:2:4:4 title "real"
#pause -1
#
#splot "0_ep.data"  u 1:2:5:5 title "imaginary"
#pause -1
#
#splot "0_ep.data"  u 1:2:6 title "absolute"
#pause -1

set pm3d
set term png
set output "1r.png"
splot "0-ep.data"  u 1:2:4:4 title "real"

set output "1r.png"
splot "0-ep.data"  u 1:2:5:5 title "imaginary"
pause -1

splot "0-ep.data"  u 1:2:6 title "absolute"
pause -1
