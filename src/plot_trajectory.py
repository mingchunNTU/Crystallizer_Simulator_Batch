import os
import sys
sys.path.insert(0, os.path.abspath('./module'))

from data import *
from figure import *
from utility import *
from crystal import *

#=============================

setting_dir=r"../examples/KNO3/"
type="T"

# option for type
# T: temperature
# S: supersaturation 
# G: growth rate
# B: nucleation rate
# tau: transform time
# C: concentration
# u0s ~ u4s: seeded moment
# u0n ~ u4n: nucleation moment
# f0: PDF for size=0

#=============================

tmp1=variable_read(setting_dir+"Result/trajectory.csv")
time=tmp1[0]
name=["T","S","G","B","tau","C","u0s","u1s","u2s","u3s","u4s","u0n","u1n","u2n","u3n","u4n","f0"]
trajectory=tmp1[name.index(type)+1]

xlim=[0,0]
ylim=[0,0]
form="o--"

auto_plot(time,trajectory,xlim,ylim,form)