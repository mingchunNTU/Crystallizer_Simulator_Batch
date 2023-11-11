from data import *
from figure import *
import numpy as np 

# ==============================================================

setting_dir=r"../examples/KNO3/"

T_upper=40 # oC
T_lower=20 # oC

# ==============================================================

file=setting_dir+"physical_property.csv"
tmp1=parameter_read(file)
A=get_variable(tmp1,"A").value[0]
B=get_variable(tmp1,"B").value[0]
C=get_variable(tmp1,"C").value[0]

T=np.linspace(T_lower,T_upper)
Csat=[]
for i in range(len(T)):
    tmp2=A+B*T[i]+C*T[i]*T[i]
    Csat.append(tmp2)

x=T
y=Csat
xlabel="Temperature ($^oC$)"
ylabel="Solubility (g/L)"
title="off"
xlim=[0,0]
ylim=[0,0]
form="o--"


plot_variable(x,y,xlabel,ylabel,title,xlim,ylim,form)