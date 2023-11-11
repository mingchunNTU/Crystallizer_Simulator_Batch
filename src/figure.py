from data import *
import matplotlib.pyplot as plt
import matplotlib.dates as mdates # setting the x axis formatter

def plot_variable(x,y,xlabel,ylabel,title,xlim,ylim,form):
	"""
	Plot the x-y figure and set the format automatically
	
	:param x: data for x-axis
	:type x: list
	:param y: data for y-axis
	:type y: list
	:param xlabel: x label of the plot
	:type xlabel: string
	:param ylabel: y label of the plot
	:type ylabel: string
	:param title: title of the plot. If title="off", there's no title
	:type title: string
	:param xlim: range for x-axis. If xlim=[0,0], the default range is used
	:type xlim: list
	:param ylim: range for y-axis. If ylim=[0,0], the default range is used
	:type ylim: list
	:param form: dot form
	:type form: char
	"""

	#plt.rcParams['font.sans-serif']=['Noto Sans CJK TC']
	#plt.rcParams['axes.unicode_minus'] = False
	plt.figure(figsize=(8,6),dpi=300)
	plt.plot(x,y,form)
	plt.xlabel(xlabel,fontsize=12)
	plt.ylabel(ylabel,fontsize=12)
	plt.xticks(fontsize=10)
	plt.yticks(fontsize=10)
	# plt.gca().xaxis.set_major_locator(plt.MultipleLocator(2))

    
	if title!="off":
		plt.title(title,fontsize=18)    
	if xlim[0]!=0 or xlim[1]!=0:
		plt.xlim(xlim[0],xlim[1])
	if ylim[0]!=0 or ylim[1]!=0:
		plt.ylim(ylim[0],ylim[1])
	
		
def plot_variable_legend(x,y,legend,xlabel,ylabel,title,xlim,ylim,form):
	"""
	Plot the x-y plot with multiple legend and set the format automatically
	
	:param x: data for x-axis
	:type x: list
	:param y: data for y-axis
	:type y: list
	:param legend: legends for each line
	:type legend: list
	:param xlabel: x label of the plot
	:type xlabel: string
	:param ylabel: y label of the plot
	:type ylabel: string
	:param title: title of the plot. If title="off", there's no title
	:type title: string
	:param xlim: range for x-axis. If xlim=[0,0], the default range is used
	:type xlim: list
	:param ylim: range for y-axis. If ylim=[0,0], the default range is used
	:type ylim: list
	:param form: dot form
	:type form: char
	"""
	
	#plt.rcParams['font.sans-serif']=['Noto Sans CJK TC']
	#plt.rcParams['axes.unicode_minus'] = False
	plt.figure(figsize=(10,6),dpi=300)
	for i in range(len(legend)):
		plt.plot(x[i],y[i],form,label=legend[i])
	plt.legend(loc="best",fontsize=14)
	plt.xlabel(xlabel,fontsize=18)
	plt.ylabel(ylabel,fontsize=18)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	# plt.gca().xaxis.set_major_locator(plt.MultipleLocator(2))

    
	if title!="off":
		plt.title(title,fontsize=18)    
	if xlim[0]!=0 or xlim[1]!=0:
		plt.xlim(xlim[0],xlim[1])
	if ylim[0]!=0 or ylim[1]!=0:
		plt.ylim(ylim[0],ylim[1])
		
def plot_validation(x,y,data_range,xlabel,ylabel,title):
	"""
	Plot the model validation graph and set the format automatically
	
	:param x: data for x-axis, usually it's experimental result
	:type x: list
	:param y: data for y-axis, usually it's simulation result
	:type y: list
	:param data_range: the range for x and y axis
	:type data_range: list
	:param xlabel: x label of the plot
	:type xlabel: string
	:param ylabel: y label of the plot
	:type ylabel: string
	:param title: title of the plot. If title="off", there's no title
	:type title: string
	"""
	#plt.rcParams['font.sans-serif']=['Noto Sans CJK TC']
	#plt.rcParams['axes.unicode_minus'] = False
	plt.plot(x,y,"o")
	plt.plot(data_range,data_range,"--")
	plt.xlim(data_range[0],data_range[1])
	plt.ylim(data_range[0],data_range[1])
	plt.xlabel(xlabel,fontsize=18)
	plt.ylabel(ylabel,fontsize=18)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	if title!="off":
		plt.title(title,fontsize=18)
		
def plot_variable_time(variable_time,ylabel,title,ylim,form,interval):
	"""
	Plot the x-t figure and set the format automatically
	
	:param variable_time: the time-dependent variable to be plotted
	:type variable_time: variable_time
	:param ylabel: y label of the plot
	:type ylabel: string
	:param title: title of the plot. If title="off", there's no title
	:type title: string
	:param ylim: range for y-axis. If ylim=[0,0], the default range is used
	:type ylim: list
	:param form: dot form
	:type form: char
	:param interval: the time step of the variable_time (in unit of hour)
	:type interval: float
	"""
	date=variable_time.date
	y=variable_time.value
	#plt.rcParams['font.sans-serif']=['Noto Sans CJK TC']
	#plt.rcParams['axes.unicode_minus'] = False 
	plt.figure(figsize=([8,8]),dpi=300)
	plt.plot(date,y,form,markersize=2)
	if ylim[0]!=0 or ylim[1]!=0:
		plt.ylim(ylim[0],ylim[1])
	plt.ylabel(ylabel,fontsize=18)
	plt.yticks(fontsize=14)
  
	interval_length=int(len(variable_time.date)*interval) 
	tick_number=5
	if interval_length<6:
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
		plt.gcf().autofmt_xdate()
	elif interval_length < 150:
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
		plt.gcf().autofmt_xdate()
		plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=int(interval_length/tick_number)))
	else:
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
		plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=int(interval_length/tick_number)))
	plt.xticks(fontsize=14)
	if title!="off":
		plt.title(title,fontsize=18)
		
def plot_variable_time_legend(variable_time_list,legend,ylabel,title,ylim,form,interval):
	"""
	Plot the x-t figure with multiple legends and set the format automatically
	
	:param variable_time: the list of the time-dependent variables to be plotted
	:type variable_time: list
	:param legend: legend for each lines
	:type legend: list
	:param ylabel: y label of the plot
	:type ylabel: string
	:param title: title of the plot. If title="off", there's no title
	:type title: string
	:param ylim: range for y-axis. If ylim=[0,0], the default range is used
	:type ylim: list
	:param form: dot form
	:type form: char
	:param interval: the time step of the variable_time (in unit of hour)
	:type interval: float
	"""
	#plt.rcParams['font.sans-serif']=['Noto Sans CJK TC']
	#plt.rcParams['axes.unicode_minus'] = False 
	plt.figure(figsize=([8,8]),dpi=300)
	for i in range(len(variable_time_list)):
		date=variable_time_list[i].date
		y=variable_time_list[i].value
		plt.plot(date,y,form,markersize=2,label=legend[i])
	if ylim[0]!=0 or ylim[1]!=0:
		plt.ylim(ylim[0],ylim[1])
	plt.legend(loc="best",fontsize=14)
	plt.ylabel(ylabel,fontsize=18)
	plt.yticks(fontsize=14)
    
	interval_length=int(len(variable_time_list[0].date)*interval) 
	tick_number=5
	if interval_length<6:
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
		plt.gcf().autofmt_xdate()
	elif interval_length < 150:
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
		plt.gcf().autofmt_xdate()
		plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=int(interval_length/tick_number)))
	else:
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
		plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=int(interval_length/tick_number)))
	plt.xticks(fontsize=14)
	if title!="off":
		plt.title(title,fontsize=18)
		
def twin_plot(variable_time1,variable_time2,ylim1,ylim2):
	"""
	Plot two time-dependent variables on different y-axis scale and set the format automatically
	
	:param variable_time1: first time-dependent variable
	:type variable_time1: variable_time
	:param variable_time2: second time-dependent variable
	:type variable_time2: variable_time
	:param ylim1: the range of y-axis limit for the first time-dependent variable
	:type ylim1: list
	:param ylim2: the range of y-axis limit for the second time-dependent variable
	:type ylim2: list
	
	"""
	# plot the result
	#plt.rcParams['font.sans-serif']=['Noto Sans CJK TC']
	#plt.rcParams['axes.unicode_minus'] = False 
	plt.figure(figsize=([8,6.5]),dpi=300)
    
	[fig, ax1]=plt.subplots(figsize=[8,6.5],dpi=300)
	color1='tab:red'
	ax1.set_ylabel(variable_time1.unit,color=color1,fontsize=16)
	p1,=ax1.plot(variable_time1.date,variable_time1.value,'o',color=color1,markersize=2,label=variable_time1.name)
	if ylim1[0]!=0 or ylim1[1]!=0:
		ax1.set_ylim(ylim1[0],ylim1[1])
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
	plt.gcf().autofmt_xdate()
	ax1.tick_params(axis='y',labelcolor=color1,labelsize=16)
	ax1.tick_params(axis='x',labelsize=12)
    
	ax2=ax1.twinx()
	color2='tab:blue'
	ax2.set_ylabel(variable_time2.unit,color=color2,fontsize=16)
	p2,=ax2.plot(variable_time2.date,variable_time2.value,'o',color=color2,markersize=2,label=variable_time2.name)
	if ylim2[0]!=0 or ylim2[1]!=0:
		ax2.set_ylim(ylim2[0],ylim2[1])
	ax2.tick_params(axis='y',labelcolor=color2,labelsize=16)
	ax1.legend(handles=[p1,p2],loc=2)    
	fig.tight_layout()    
