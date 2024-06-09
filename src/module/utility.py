import numpy as np
from scipy.integrate import quad, odeint
from data import *
from figure import *

def transform(x_list,y_list,x_coordinate):
	"""
	Transform the discretized PDF to continuous PDF, which makes it easier to integrate
	
	:param x_list: crystal size (um)
	:type x_list: list
	:param y_list: density function (1/um^4)
	:type y_list: list
	:param x_coordinate: x coordinate (um)
	:type x_coordinate: float
	:return: y coordinate (1/um^4)
	:rtype: float
	"""
	
	y_coordinate=0.0
	if x_coordinate < x_list[0]:
		y_coordinate=y_list[0]
	elif x_coordinate > x_list[-1]:
		y_coordinate=y_list[-1]
	for i in range(len(x_list)-1):
		if x_coordinate == x_list[i]:
			y_coordinate=y_list[i]
		elif x_coordinate > x_list[i] and x_coordinate < x_list[i+1]:
			y_coordinate=y_list[i]*(x_list[i+1]-x_coordinate)/(x_list[i+1]-x_list[i])+y_list[i+1]*(x_coordinate-x_list[i])/(x_list[i+1]-x_list[i])
		elif x_coordinate == x_list[i+1]:
			y_coordinate=y_list[i+1]
	
	return y_coordinate

def moment_calculation(size,density_function,order):
	"""
	Calculate the desired order moment of the specified CSD
	
	:param size: crystal size (um)
	:type size: list
	:param density_function: density function (1/um^4)
	:type density_function: list
	:param order: moment order
	:type order: float
	:return: moment
	:rtype: float 
	"""
	
	x_list=size
	y_list=density_function
	continuous_function=lambda x: transform(x_list,y_list,x)*x**order
	output=quad(continuous_function,size[0],size[-1])[0]
	
	return output

def list_operation(input_list,operation,value):
	"""
	Shift or Scale the CSD by a specified value
	
	:param input_list: input list
	:type input_list: list
	:param operation: + for shift, * for scale
	:type operation: char
	:param value: shift or scale constant
	:type value: float
	:return: processed list
	:rtype: list
	"""
	
	output=[]
	for i in range(len(input_list)):
		if operation == "+":
			output.append(input_list[i]+value)
		elif operation == "*":
			output.append(input_list[i]*value)
	return output

def CSD_normalization(solid_loading,crystal_density,shape_factor,size,crude_density_function):
	"""
	Convert the crude CSD to a CSD that fits the specified solid loading
	
	:param solid_loading: solid loading of the crystal (g/L)
	:type solid_loading: float
	:param crystal_density: crystal density (g/L)
	:type crystal_density: float
	:param shape_factor: shape factor (-)
	:type shape_factor: float
	:param size: crystal size (um)
	:type size: list
	:param crude_density_function: crude density function (-)
	:type crude_density_function: list
	:return: return the size (list) and density function (list)
	"""
	
	
	normalization_constant=solid_loading/crystal_density/shape_factor/moment_calculation(size,crude_density_function,3)
	density_function=list_operation(crude_density_function,'*',normalization_constant)
	density_function2=variable("density function","$1/\mu m^4$",density_function)
	
	return size, density_function2
	
def generate_seed_parabolic(seed_loading,crystal_density,shape_factor,mesh_size,r0,w):
	"""
	Generate the seed CSD (size and density function) through seed usage
	
	:param seed_loading: seed loading (g/L)
	:type seed_loading: float
	:param crystal_density: crystal density (g/L)
	:type crystal_density: float
	:param shape_factor: shape factor (-)
	:type shape_factor: float
	:param mesh_size: mesh size (um)
	:type mesh_size: float
	:param r0: seed mean size (um)
	:type r0: float
	:param w: seed CSD width (-)
	:type w: float
	:return: size and density function
	:rtype: list
	"""
	
	mesh_number=int(2*r0*w/mesh_size)+1
	size=np.linspace(r0*(1-w),r0*(1+w),mesh_number)
	parabolic=lambda x: (x-r0*(1-w))*(r0*(1+w)-x)
	crude_density_function=[]
	for i in range(len(size)):
		crude_density_function.append(parabolic(size[i]))
	
	CSD_normalization(seed_loading,crystal_density,shape_factor,size,crude_density_function)
	size, density_function=CSD_normalization(seed_loading,crystal_density,shape_factor,size,crude_density_function)
	
	return size, density_function	

def generate_temperature_trajectory(batch_time,time_step,initial_T,end_T,order):
	"""
	Generate the temperature trajectory with specified order
	
	:param batch_time: batch time (min)
	:type batch_time: float
	:param time_step: time step (min)
	:type time_step: float
	:param initial_T: initial temperature (C)
	:type initial_T: float
	:param end_T: end temperature (C
	:type end_T: float
	:param order: order of cooling trajectory (-)
	:type order: float
	:return: time and temperature
	:rtype: list
	
	"""
	
	mesh_number=int(batch_time/time_step)+1
	time=np.linspace(0,batch_time,mesh_number)
	temperature=[]
	for i in range(mesh_number):
		tmp1=initial_T-(initial_T-end_T)*(time[i]/batch_time)**order
		temperature.append(tmp1)
	
	time2=variable("time","min",time)
	temperature2=variable("temperature","$^oC$",temperature)


	return time2, temperature2

def DF_to_NF(size,density_function,mesh_size):
	"""
	Transform the number density function to number fraction
	
	:param size: size
	:type size: variable
	:param density_function: number density function
	:type density_function: variable
	:param mesh_size: mesh size of the size distribution (um)
	:type mesh_size: float
	:return: size and number fraction
	:rtype: variable
	"""
	
	size_output=variable("size","$\mu m$",[])
	number_fraction=variable("number fraction","-",[])
	
	mesh_number=int((size.value[-1]-size.value[0])/mesh_size)+1
	tmp1=np.linspace(size.value[0],size.value[-1],mesh_number)
	normalization_constant=moment_calculation(size.value,density_function.value,0)
	
	for i in range(len(tmp1)-1):
		L_mean=(tmp1[i]+tmp1[i+1])/2
		f1=transform(size.value,density_function.value,tmp1[i])
		f2=transform(size.value,density_function.value,tmp1[i+1])
		tmp2=(f1+f2)/2*mesh_size/normalization_constant
		size_output.value.append(L_mean)
		number_fraction.value.append(tmp2)
		
	return size_output, number_fraction

def NF_to_VF(size,number_fraction):
	"""
	Transform the number fraction to volume fraction
	
	:param size: size
	:type size: variable
	:param number_fraction: number fraction
	:type number_fraction: variable
	:return: size and volume fraction
	:rtype: variable
	"""
	
	volume_fraction=variable("volume fraction","-",[])
	
	total_volume=0
	for i in range(len(size.value)):
		tmp1=size.value[i]**3*number_fraction.value[i]
		total_volume=total_volume+tmp1
	for i in range(len(size.value)):
		tmp1=size.value[i]**3*number_fraction.value[i]
		volume_fraction.value.append(tmp1/total_volume)	
	
	return size,volume_fraction

def VF_to_NF(size,volume_fraction):
	"""
	Transform the number fraction to volume fraction
	
	:param size: size
	:type size: variable
	:param volume_fraction: volume fraction
	:type volume_fraction: variable
	:return: size and number fraction
	:rtype: variable
	"""
	
	number_fraction=variable("number fraction","-",[])
	
	total_number=0
	for i in range(len(size.value)):
		tmp1=volume_fraction.value[i]/size.value[i]**3
		total_number=total_number+tmp1
	for i in range(len(size.value)):
		tmp1=volume_fraction.value[i]/size.value[i]**3
		number_fraction.value.append(tmp1/total_number)	
	
	return size,number_fraction

def auto_plot(variable_x,variable_y,xlim,ylim,form):
	"""
	Plot the variables and generate xlabel and ylabel automatically
	
	:param variable_x: data for x-axis
	:type variable_x: variable
	:param variable_y: data for y-axis
	:type variable_y: variable
	:param title: title of the plot. If title="off", there's no title
	:type title: string
	:param xlim: range for x-axis. If xlim=[0,0], the default range is used
	:type xlim: list
	:param ylim: range for y-axis. If ylim=[0,0], the default range is used
	:type ylim: list
	:param form: dot form
	:type form: char
	"""
	
	x=variable_x.value
	y=variable_y.value
	xlabel=variable_x.name+" ("+variable_x.unit+")"
	ylabel=variable_y.name+" ("+variable_y.unit+")"
	title="off"
	
	plot_variable(x,y,xlabel,ylabel,title,xlim,ylim,form)
	
