from data import *
from figure import *
from utility import *
from crystal import *
import sys


tmp1=crystal()
tmp1.initialize("physical_property.csv")

# read the operating condition
operating_condition=parameter_read("operating_condition.csv")
physical_property=parameter_read("physical_property.csv")

# set the seed condition
seed_loading=get_variable(operating_condition,"seed_loading")
r0=get_variable(operating_condition,"r0")
w=get_variable(operating_condition,"w")
mesh_size=get_variable(operating_condition,"mesh_size")
crystal_density=get_variable(physical_property,"density")
shape_factor=get_variable(physical_property,"shape_factor")


size,density_function=generate_seed_parabolic(seed_loading.value[0],crystal_density.value[0],shape_factor.value[0],mesh_size.value[0],r0.value[0],w.value[0])
tmp1.set_seed(size,density_function)

# set the temperature trajectory
batch_time= get_variable(operating_condition,"batch_time")
time_step=get_variable(operating_condition,"time_step")
initial_T=get_variable(operating_condition,"initial_T")
end_T=get_variable(operating_condition,"end_T")
order=get_variable(operating_condition,"order")

time, temperature=generate_temperature_trajectory(batch_time.value[0],time_step.value[0],initial_T.value[0],end_T.value[0],1)
tmp1.set_trajectory(time, temperature)
	
# solve CSD and trajectory
tmp1.solve()

# plot and save the result

if sys.argv[1]=='--get_CSD':
	size, density_function=tmp1.get_CSD()
	size, number_fraction=DF_to_NF(size,density_function,float(sys.argv[2]))
	
	if sys.argv[3]=='number':
		plot_variable_II(size,number_fraction,"off",[0,0],[0,0],'o')
		plt.show()
	
		output=[size,number_fraction]
		variable_output(output,"CSD.csv")
		
	if sys.argv[3]=='volume':
		size,volume_fraction=NF_to_VF(size,number_fraction)
		plot_variable_II(size,volume_fraction,"off",[0,0],[0,0],'o')
		plt.show()
	
		output=[size,volume_fraction]
		variable_output(output,"CSD.csv")
	
elif sys.argv[1]=='--get_trajectory':
	time,trajectory=tmp1.get_trajectory(sys.argv[2])
	plot_variable_II(time,trajectory,"off",[0,0],[0,0],'-')
	plt.show()
	
	output=[time,trajectory]
	variable_output(output,"trajectory.csv")


