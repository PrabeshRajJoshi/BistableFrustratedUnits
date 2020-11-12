'''
Python file to carry out the simulation
'''

import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Process
from subprocess import call
import datetime
import random

# make time stamp for file handling
now = datetime.date.now()
now_date = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
now_time = str(now.hour).zfill(2) + str(now.minute).zfill(2) + str(now.second).zfill(2)

# options to execute the fortran script
nx = "004"
bR = "0.0"
rt = "1000.0"
dt = "0.01"
a0 = "18.0"
a1 = "0.0" # insert non-zero value (15.0) if increasing alpha is needed
sa = "4000.0" # speed of alpha
ws = "1.0"
sA = "1" # seeds
sB = "2"

file_prefix = nx + "x" + nx + "_" + now_date + now_time
fn = file_prefix

#fortran file to execute (without extension!)
fortran_script = "cbfu"

#name of the directory where movie and data files hould transfer to
directory_name = file_prefix + "_alpha0_"+a0

#execute the script where all local functions are defined
execfile("cbfu_modules.py")

# call function to execute fortran script
execute_fortran(fortran_script, nx, bR, rt, dt, a0, a1, sa, ws, sA, sB, fn)

# call function to make images
#   file_length as the third argument to avoid memory error
#   "A"/"B" as the fourth argument to plot for data of A/B
#   "yes"/"no" as the fifth argument to make/not make a movie.

# writing step to make images accordingly
save_step = int(float(ws))
# time to stop simulation video
stop = 1000
imaging_and_movie(file_prefix, save_step, "A", "yes", stop)

# call function to make plots of segments of time
plt_title = "$ \\ beat_R: "+ bR + "\; \\alpha_0: "+a0 + "\; \Delta\\alpha: "+a1 + "\; \\alpha_{speed}: "+sa + "\; configuration: "+fortran_script + "$"
nbfu = int(nx) * int(nx)

# define the segment where the plot is needed
start = 0
stop = 1000
segment_plots(file_prefix, "A", nbfux, ws, start, stop, plt_title)

# call function to organize files (directory_name) is moved to either "changing_alpha" or "constant_alpha")
organize_files(directory_name, file_prefix, fortran_script, a1, bR)

print("Task finished!\n")