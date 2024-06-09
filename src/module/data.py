import numpy as np
import csv
import copy
from datetime import *

def csv_input(file):
	"""
	Read the csv file as a list
    
	:param file: csv file name
	:type file: string
	:return: content of csv file
	:rtype: list
    
	"""
    
	tmp=open(file,'r',newline='')
	read=csv.reader(tmp)
	output=[]
	for x in read:
		output.append(x)
	tmp.close()
	return output
	
def csv_output(target,file):
	"""
	Export the list as a csv_file
	
	:param target: the list to be exported
	:type target: list
	:param file: the csv file name
	:type file: string
	
	"""
	
	tmp=open(file,'w',newline='')
	write=csv.writer(tmp)
	for x in target:
		write.writerow(x)
	tmp.close()
	
def string_float_transform(string):
	"""
	Transform the string to float number. If the string can not transform to number, output 0 and set bad_data_checker to 1
	
	:param string: the string to be transformed
	:type string: string
	:return: The float and bad_data_checker
	:rtype: float and int 
	"""
	
	bad_data_checker=0
	if string[0].isdigit() or string[0]=="-":
		output=float(string)
	else:
		output=float(0)
		bad_data_checker=1
	return output,bad_data_checker
	
class variable:
	"""
	A class built for the variables that contain name and physical unit. It's defined through init function.
	
	"""	
	def __init__(self,name,unit,value):
		self.name=name
		self.unit=unit
		self.value=value
	
def variable_read(file):
	"""
	Read the csv file and transform it into list of variables. Print the warning message if the transformation is not successful
 
	:param file: csv file
	:type file: string
	
	"""
	tmp1=csv_input(file)
	output=[]
	bad_data_checker=0
	for i in range(len(tmp1[0])):
		name=tmp1[0][i]
		unit=tmp1[1][i]
		value=[]
		for j in range(len(tmp1)-2):
			tmp2,tmp3=string_float_transform(tmp1[j+2][i])
			value.append(tmp2)
			if tmp3 != 0:
				bad_data_checker=1        
		output.append(variable(name,unit,value))
	if bad_data_checker !=0:
		print("Bad data detected(not float)")
	return output
	
def parameter_read(file):
	"""
	Read the csv file that contains multiple parameters and transform it into list of variables
	
	:param file: csv file
	:type file: string
	:return: list of parameters
	:rtype: list
	"""
	
	tmp1=csv_input(file)
	iteration=int((len(tmp1)+1)/5)
	output=[]
	for i in range(iteration):
		for j in range(len(tmp1[5*i+1])):
			name=tmp1[5*i+1][j]
			unit=tmp1[5*i+2][j]
			value=[]
			if tmp1[5*i+3][j] != "":
				value.append(float(tmp1[5*i+3][j]))
			tmp2=variable(name,unit,value)
			output.append(tmp2)
	return output
	
	
def get_variable(variable_list,variable_name):
	"""
	Return the variable that fits the desired variable name from a variable list
	
	:param variable_list: variable list
	:type variable_list: list
	:param variable_name: desired variable
	:type variable_name: string
	:return: desired variable
	:rtype: variable 
	"""
	
	name_list=[]
	for i in range(len(variable_list)):
		name_list.append(variable_list[i].name)
	
	return variable_list[name_list.index(variable_name)]

	
def variable_output(variable_list,file):
	"""
	Export the variables to csv file
    
	:param variable_list: variables to be exported
	:type variable_list: list
	:param file: csv file
	:type file: string
    
	"""
	name_list=[]
	unit_list=[]
	value_list=[]
	for i in range(len(variable_list)):
		name_list.append(variable_list[i].name)
		unit_list.append(variable_list[i].unit)
		value_list.append(variable_list[i].value)
    
	output=[]
	output.append(name_list)
	output.append(unit_list)
	for i in range(len(value_list[0])):
		tmp1=[]
		for j in range(len(value_list)):
			tmp1.append(value_list[j][i])
		output.append(tmp1)

	csv_output(output,file)	
	
class variable_time:
	"""
	A class for variables that is time-dependent. It's defined through init function.
	"""
    
	def __init__(self,name,unit,value,date):
		self.name=name
		self.unit=unit
		self.value=value
		self.date=date
		self.init_time=date[0]
		self.end_time=date[-1]	
		
	def set_time(self,time1,start_hour,time2,end_hour):
        
		"""
		Set the time interval for the variable_time.
    	
		:param time1: start date in the format of %Y/%m/%d
		:type time1: string
		:param start_hour: start hour for the start date
		:type start_hour: int
		:param time2: end date in the format of %Y/%m/%d
		:type time2: string
		:param end_hour: end hour for the start date
		:type end_hour: int
    	
		"""
        
		init_time=datetime.strptime(time1,"%Y/%m/%d")
		end_time=datetime.strptime(time2,"%Y/%m/%d")
    
		index1=0
		index2=0
		for x in range(len(self.date)):
			tmp1=self.date[x]
			if tmp1.year==init_time.year and tmp1.month==init_time.month and tmp1.day==init_time.day and tmp1.hour==start_hour:
				index1=x
				break
		for x in range(len(self.date)):
			tmp1=self.date[x]
			if tmp1.year==end_time.year and tmp1.month==end_time.month and tmp1.day==end_time.day and tmp1.hour==end_hour:
				index2=x
				break    
    
		tmp_date=self.date
		tmp_value=self.value
		self.date=[]
		self.value=[]
		for x in range(index1,index2+1):
			self.date.append(tmp_date[x])
			self.value.append(tmp_value[x])
		self.init_time=init_time
		self.end_time=end_time
	
	def summary(self):
		"""
		Export the mean, standard deviation and relative standard deviation of the variable_time
    	
		:return: mean, std and rstd
		"""

		mean=np.mean(self.value)
		std=np.std(self.value)
		rstd=std/mean*100
		return mean,std,rstd
	
def variable_time_read(file,date_form):
	"""
	Read the csv file and transform it into variable_time. Print the warning message if the transformation is not successful
	
	:param file: csv file
	:type file: string
	:date_form: the date format of the csv data
	:type form: string
	"""

	tmp1=csv_input(file)
	date=[]
	for i in range(len(tmp1)-2):
		date.append(datetime.strptime(tmp1[i+2][0],date_form))
	output=[]
	bad_data_checker=0
	for i in range(len(tmp1[0])-1):
		name=tmp1[0][i+1]
		unit=tmp1[1][i+1]
		value=[]
		for j in range(len(tmp1)-2):
			tmp2,tmp3=string_float_transform(tmp1[j+2][i+1])
			value.append(tmp2)
			if tmp3 != 0:
				bad_data_checker=1
		output.append(variable_time(name,unit,value,date))
	if bad_data_checker==1:
		print("Bad data detected(not float)")
	return output	
	
	
	
