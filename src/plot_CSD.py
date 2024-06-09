import os
import sys
sys.path.insert(0, os.path.abspath('./module'))

from data import *
from figure import *
from utility import *
from crystal import *

#=============================

setting_dir=r"../examples/KNO3/"
type="number" # number or volume

#=============================

if type=="number":
    file=setting_dir+"Result/CSD(number).csv"
elif type=="volume":
    file=setting_dir+"Result/CSD(volume).csv"

tmp1=variable_read(file)
x=tmp1[0].value
y=tmp1[1].value
xlabel="size ($\mu m$)"
ylabel="fraction (-)"
title="Crystal Size Distribution ("+type+")"
xlim=[0,0]
ylim=[0,0]
form="o--"

plot_variable(x,y,xlabel,ylabel,title,xlim,ylim,form)