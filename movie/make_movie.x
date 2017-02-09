#mencoder "mf://*.png" -vf scale=800:600 -mf fps=6 -o brownian.avi -ovc lavc -lavcopts vcodec=mjpeg:vbitrate=800
#mencoder "mf://*.png" -mf fps=24 -o output.avi -ovc lavc vcodec=mjpeg:vbitrate=800
#mencoder "mf://*.png" -mf fps=6 -o output.avi -ovc lavc 
#mencoder "mf://plot*.png" -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi

mencoder mf://plot*.png -mf w=800:h=600:fps=25:type=png -ovc copy -oac copy -o output.avi

