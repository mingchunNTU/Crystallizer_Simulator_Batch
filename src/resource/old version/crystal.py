from data import *
from figure import *
from utility import *
import numpy as np
from scipy.integrate import quad, odeint
import time # calculate the execution time

class crystal():
	"""
	Module used for crystallization simulation
	"""

	def initialize(self,filename):
		"""
		Read in the physical property of the crystal
		
		:param filename: physical property file name
		:type filename: string
		""" 
		
		# read the physical property
		tmp1=parameter_read(filename)
		self.density=get_variable(tmp1,"density")
		self.shape_factor=get_variable(tmp1,"shape_factor")
		self.kb=get_variable(tmp1,"kb")
		self.b=get_variable(tmp1,"b")
		self.kg=get_variable(tmp1,"kg")
		self.g=get_variable(tmp1,"g")
		self.A=get_variable(tmp1,"A")
		self.B=get_variable(tmp1,"B")
		self.C=get_variable(tmp1,"C")
		
		# initialize the trajectory
		self.trajectory_time=variable("time","min",[])
		self.trajectory_temperature=variable("temperature","$^oC$",[])
		self.trajectory_supersaturation=variable("supersaturation","-",[])
		self.trajectory_growth=variable("growth rate","$1/(\mu m^3*min)$",[])
		self.trajectory_nucleation=variable("nucleation rate","$\mu m/s$",[])
		self.trajectory_tau=variable("transform time","$\mu m$",[]) # used for CSD calculation
		self.trajectory_f0=variable("fn_0","$1/\mu m^4$",[]) # used for CSD calculation
		
		# initialize CSD
		self.size=variable("size","$\mu m$",[])
		self.density_function=variable("density function","$1/\mu m^4$",[])
		
		
		
	def solubility(self,T):
		"""
		Calculate the saturated concentration under specified temperature
		
		:param T: temperature (C)
		:type T: float
		:return: saturated concentration (g/L)
		:rtype: float
		"""
		Csat=self.A.value[0]+self.B.value[0]*T+self.C.value[0]*T**2
		
		return Csat
		
	def set_trajectory(self,time,temperature):
		"""
		Set the cooling trajectory
		
		:param time: time (min)
		:type time: list
		:param temperature: temperature (C)
		:type temperature: list
		"""
		self.trajectory_time.value=time
		self.trajectory_temperature.value=temperature
		
	def set_seed(self,size,density_function):
		"""
		Set the CSD of seed
		
		:param size: crystal size (um)
		:type size: list
		:param density_function: density function (1/um^4)
		:type density_function: list
		""" 
		self.size.value=size
		self.density_function.value=density_function
		
	def moment_equation(self,y,t):
		"""
		Define the moment equation to be solved. y[0]~y[4]: seed moment (order=index). y[5]~y[9]: nucleate moment (order=index-5). y[10]: concentration 
		
		:param y: moment
		:type y: list
		:param t: time
		:type t: float
		:return: dy/dt
		:rtype: list
		"""
	
		T=transform(self.trajectory_time.value,self.trajectory_temperature.value,t)
		Csat=self.solubility(T)
		S=(y[10]-Csat)/Csat
		if S<0:
			S=0
		G=self.kg.value[0]*S**self.g.value[0]
		B=self.kb.value[0]*S**self.b.value[0]*(y[3]+y[8])
		
		dydt=[]
		dydt.append(0) # dy[0]/dt
		dydt.append(y[0]*G) # dy[1]/dt
		dydt.append(2*y[1]*G) # dy[2]/dt
		dydt.append(3*y[2]*G) # dy[3]/dt
		dydt.append(4*y[3]*G) # dy[4]/dt
		dydt.append(B) # dy[5]/dt
		dydt.append(y[5]*G) # dy[6]/dt
		dydt.append(2*y[6]*G) # dy[7]/dt
		dydt.append(3*y[7]*G) # dy[8]/dt
		dydt.append(4*y[8]*G) # dy[9]/dt
		dydt.append(-3*(y[2]+y[7])*G*self.shape_factor.value[0]*self.density.value[0])
		
		return dydt
		
	def solve(self):
		"""
		Solve the moment equation, get various trajectory and final CSD of the crystallizer
		
		"""
		
		start_time=time.time() # keep track of the execution time
		
		# calculate the initial condition for moment equation
		initial=[]
		for i in range(5):
			initial.append(moment_calculation(self.size.value,self.density_function.value,i))
			
		for i in range(5):
			initial.append(0.0)
			
		initial_T=self.trajectory_temperature.value[0]
		C0=self.solubility(initial_T)
		initial.append(C0)
		
		# solve the moment equation
		self.trajectory_moment=odeint(self.moment_equation,initial,self.trajectory_time.value)
		
		# update the trajectory of moment
		self.trajectory_u0s=variable("$\mu_{s,0}$","$1/\mu m^3$",self.trajectory_moment[:,0])
		self.trajectory_u1s=variable("$\mu_{s,1}$","$1/\mu m^2$",self.trajectory_moment[:,1])
		self.trajectory_u2s=variable("$\mu_{s,2}$","$1/\mu m$",self.trajectory_moment[:,2])
		self.trajectory_u3s=variable("$\mu_{s,3}$","-",self.trajectory_moment[:,3])
		self.trajectory_u4s=variable("$\mu_{s,4}$","$\mu m$",self.trajectory_moment[:,4])
		self.trajectory_u0n=variable("$\mu_{n,0}$","$1/\mu m^3$",self.trajectory_moment[:,5])
		self.trajectory_u1n=variable("$\mu_{n,1}$","$1/\mu m^2$",self.trajectory_moment[:,6])
		self.trajectory_u2n=variable("$\mu_{n,2}$","$1/\mu m$",self.trajectory_moment[:,7])
		self.trajectory_u3n=variable("$\mu_{n,3}$","-",self.trajectory_moment[:,8])
		self.trajectory_u4n=variable("$\mu_{n,4}$","$\mu m$",self.trajectory_moment[:,9])
		self.trajectory_C=variable("concentration","g/L",self.trajectory_moment[:,10])
		
		# update the S,G,B trajectory
		C=self.trajectory_moment[:,-1]
		u3s=self.trajectory_moment[:,3]
		u3n=self.trajectory_moment[:,8]
		for i in range(len(C)):
			T=self.trajectory_temperature.value[i]
			Csat=self.solubility(T)
			S=(C[i]-Csat)/Csat
			if S<0:
				S=0
			G=self.kg.value[0]*S**self.g.value[0]
			B=self.kb.value[0]*S**self.b.value[0]*(u3s[i]+u3n[i])
			
			self.trajectory_supersaturation.value.append(S)
			self.trajectory_growth.value.append(G)
			self.trajectory_nucleation.value.append(B)
			
		# update tau and f0 trajectory for CSD calculation
		self.trajectory_tau.value.append(0)
		self.trajectory_f0.value.append(0)
		for i in range(len(self.trajectory_time.value)-1):
			delta_t=self.trajectory_time.value[i+1]-self.trajectory_time.value[i]
			G_mean=(self.trajectory_growth.value[i+1]+self.trajectory_growth.value[i+1])/2
			tau_old=self.trajectory_tau.value[i]
			G_ins=self.trajectory_growth.value[i+1]
			B_ins=self.trajectory_nucleation.value[i+1]
			self.trajectory_tau.value.append(tau_old+G_mean*delta_t)
			if G_ins==0:
				self.trajectory_f0.value.append(0)
			else:
				self.trajectory_f0.value.append(B_ins/G_ins)
		
		# calculate the CSD
		size_end=[]
		density_function_end=[]
		tau_end=self.trajectory_tau.value[-1]
		
		# calculate the CSD of nucleated crystal
		for i in range(len(self.trajectory_tau.value)):
			index=len(self.trajectory_tau.value)-i-1
			size_end.append(tau_end-self.trajectory_tau.value[index])
			density_function_end.append(self.trajectory_f0.value[index])
			
		# calculate the CSD of seeded crystal
		size_seed=list_operation(self.size.value,"+",tau_end)
		size_end.extend(size_seed)
		density_function_end.extend(self.density_function.value)
		
		self.size_end=variable("size","$\mu m$",size_end)
		self.density_function_end=variable("density function","$1/\mu m^4$",density_function_end)
			
		# print the calculation time
		end_time=time.time()
		program_time=round(end_time-start_time,0)
		string="execution time= "+str(program_time)+" s"
		print(string)
		
	def get_CSD(self):
		"""
		Return the size and density function of the final crystal
		
		:return: size (variable) and density function (variable)
		:rtype: variable
		"""
		
		return self.size_end, self.density_function_end
		
	def get_trajectory(self,name):
		"""
		Return the desired trajectory
		
		:param name: desired trajectory name (T,S,G,B,tau,C,u0s~u4s,u0n~u4n)
		:type name: string
		:return: time (variable) and trajectory (variable)
		:rtype: variable
		"""
		
		if name=='T':
			return self.trajectory_time,self.trajectory_temperature
		elif name=='S':
			return self.trajectory_time,self.trajectory_supersaturation
		elif name=='G':
			return self.trajectory_time,self.trajectory_growth
		elif name=='B':
			return self.trajectory_time,self.trajectory_nucleation
		elif name=='tau':
			return self.trajectory_time,self.trajectory_tau
		elif name=='C':
			return self.trajectory_time,self.trajectory_C
		elif name=='u0s':
			return self.trajectory_time,self.trajectory_u0s
		elif name=='u1s':
			return self.trajectory_time,self.trajectory_u1s
		elif name=='u2s':
			return self.trajectory_time,self.trajectory_u2s
		elif name=='u3s':
			return self.trajectory_time,self.trajectory_u3s
		elif name=='u4s':
			return self.trajectory_time,self.trajectory_u4s
		elif name=='u0n':
			return self.trajectory_time,self.trajectory_u0n
		elif name=='u1n':
			return self.trajectory_time,self.trajectory_u1n
		elif name=='u2n':
			return self.trajectory_time,self.trajectory_u2n
		elif name=='u3n':
			return self.trajectory_time,self.trajectory_u3n
		elif name=='u4n':
			return self.trajectory_time,self.trajectory_u4n
		
		
		
		
		
		
		
			
			
		
	
