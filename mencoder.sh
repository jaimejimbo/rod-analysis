#!/bin/bash

# Creates a video of all images JPG.

mencoder mf://*.JPG -mf fps=15:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o name.avi
