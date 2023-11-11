from data import *
from figure import *
from utility import *
from crystal import *

#=============================

setting_dir=r"../examples/KNO3/"
CSD_mesh_size=1

#=============================

# read the crystal property
tmp1=crystal()
tmp1.read_property(setting_dir+"physical_property.csv")


# set the operating condition
operating_condition=parameter_read(setting_dir+"operating_condition.csv")
physical_property=parameter_read(setting_dir+"physical_property.csv")

# seed condition
seed_loading=get_variable(operating_condition,"seed_loading")
r0=get_variable(operating_condition,"r0")
w=get_variable(operating_condition,"w")
mesh_size=get_variable(operating_condition,"mesh_size")
crystal_density=get_variable(physical_property,"density")
shape_factor=get_variable(physical_property,"shape_factor")


size,density_function=generate_seed_parabolic(seed_loading.value[0],crystal_density.value[0],shape_factor.value[0],mesh_size.value[0],r0.value[0],w.value[0])
tmp1.set_seed(size,density_function)

# temperature trajectory
batch_time= get_variable(operating_condition,"batch_time")
time_step=get_variable(operating_condition,"time_step")
initial_T=get_variable(operating_condition,"initial_T")
end_T=get_variable(operating_condition,"end_T")
order=get_variable(operating_condition,"order")

time, temperature=generate_temperature_trajectory(batch_time.value[0],time_step.value[0],initial_T.value[0],end_T.value[0],order.value[0])
tmp1.set_cooling_trajectory(time, temperature)

# solve and output CSD and trajectory
tmp1.solve()
tmp1.output(setting_dir,CSD_mesh_size)
