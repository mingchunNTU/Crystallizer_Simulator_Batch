from data import *
from figure import *
from utility import *
import numpy as np
from scipy.integrate import quad, odeint
import time as tm # calculate the execution time

class CSD:
    """
    Module used to represent crystal size distribution
    """

    def __init__(self,size,fraction):
        self.size=size
        self.fraction=fraction

class PDF:

    """
    Module used to represent population density function
    """
    def __init__(self,size,density_function):
        self.size=size
        self.density_function=density_function

class trajectory:
    """
    Module used to represent trajectory

    """

    def __init__(self,time,trajectory_in):
        self.time=time
        self.trajectory=trajectory_in

class crystal:
    """
    Module used for crystallization simulation
    """

    def __init__(self):
        
        # initialize the trajectory 
        time=variable("time","min",[])
        T=variable("temperature","$^oC$",[])
        S=variable("supersaturation","-",[])
        G=variable("growth rate","$\mu m/min$",[])
        B=variable("nucleation rate","$1/(\mu m^3*min)$",[])
        tau=variable("transform time","\mu m",[])
        C=variable("concentration","g/L",[])
        u0s=variable("$\mu_{s,0}$","$1/\mu m^3$",[])
        u1s=variable("$\mu_{s,1}$","$1/\mu m^2$",[])
        u2s=variable("$\mu_{s,2}$","$1/\mu m$",[])
        u3s=variable("$\mu_{s,3}$","-",[])
        u4s=variable("$\mu_{s,4}$","$\mu m$",[])
        u0n=variable("$\mu_{n,0}$","$1/\mu m^3$",[])
        u1n=variable("$\mu_{n,1}$","$1/\mu m^2$",[])
        u2n=variable("$\mu_{n,2}$","$1/\mu m$",[])
        u3n=variable("$\mu_{n,3}$","-",[])
        u4n=variable("$\mu_{n,4}$","$\mu m$",[])
        f0=variable("PDF for size=0","$1/\mu m^4$",[])
        
        tmp1=[T,S,G,B,tau,C,u0s,u1s,u2s,u3s,u4s,u0n,u1n,u2n,u3n,u4n,f0]
        self.trajectory_name=["T","S","G","B","tau","C","u0s","u1s","u2s","u3s","u4s","u0n","u1n","u2n","u3n","u4n","f0"]
        self.trajectory_list=[]
        for i in range(len(tmp1)):
            self.trajectory_list.append(trajectory(time,tmp1[i]))

        # initialize the PDF and CSD
        size=variable("size","$\mu m$",[])
        density_function=variable("density function","$1/\mu m^4$",[])
        self.seed_PDF=PDF(size,density_function)
        self.product_PDF=PDF(size,density_function)
    
    def get_trajectory_index(self,name):
        
        """
        Get the trajectory index of specify name

        :param name: trajectory name
        :type name: string
        :return: list index
        :rtype: int
        """
        index="none"
        for i in range(len(self.trajectory_list)):
            if self.trajectory_list[i].trajectory.name==name:
                index=i
        if index=="none":
            print("the trajectory name is wrong")
        else:
            return index
        
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

    def read_property(self,filename):
        """
        Read the physical property of the crystal
    
        :param filename: physical property file name
        :type filename: string
        """
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
    
    def set_cooling_trajectory(self,time,temperature):
        """
        Set the cooling trajectory
        
        :param time: time (min)
        :type time: list
        :param temperature: temperature (C)
        :type temperature: list
        """
        self.trajectory_list[self.get_trajectory_index("temperature")].time=time
        self.trajectory_list[self.get_trajectory_index("temperature")].trajectory=temperature

    def set_seed(self,size,density_function):
        """
        Set the PDF of seed
        
        :param size: crystal size (um)
        :type size: list
        :param density_function: density function (1/um^4)
        :type density_function: list
        """ 
        self.seed_PDF.size.value=size
        self.seed_PDF.density_function=density_function

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

        # calculate supersaturation, nucleation rate and growth rate
        time_list=self.trajectory_list[self.get_trajectory_index("temperature")].time.value
        temperature_list=self.trajectory_list[self.get_trajectory_index("temperature")].trajectory.value
        T=transform(time_list,temperature_list,t)
        Csat=self.solubility(T)
        S=(y[10]-Csat)/Csat
        if S<0:
            S=0
        G=self.kg.value[0]*S**self.g.value[0]
        B=self.kb.value[0]*S**self.b.value[0]*(y[3]+y[8])
        
        # define the moment_equation
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
        Solve the moment equation, get various trajectory and final PDF of the crystallizer
        
        """
        
        start_time=tm.time() # keep track of the execution time
        
        # calculate the initial condition for moment equation
        initial=[]
        for i in range(5):
            size=self.seed_PDF.size.value
            density_function=self.seed_PDF.density_function.value
            initial.append(moment_calculation(size,density_function,i))
            
        for i in range(5):
            initial.append(0.0)
        
        time=self.trajectory_list[self.get_trajectory_index("temperature")].time.value
        T=self.trajectory_list[self.get_trajectory_index("temperature")].trajectory.value

        initial_T=T[0]
        C0=self.solubility(initial_T)
        initial.append(C0)
        
        # solve the moment equation
        moment_trajectory=odeint(self.moment_equation,initial,time)
        
        # update the trajectory of moment
        for i in range(10):
            self.trajectory_list[i+6].time.value=time
            self.trajectory_list[i+6].trajectory.value=moment_trajectory[:,i]

        self.trajectory_list[self.get_trajectory_index("concentration")].time.value=time
        self.trajectory_list[self.get_trajectory_index("concentration")].trajectory.value=moment_trajectory[:,10]
        
        # update the S,G,B trajectory
        self.trajectory_list[self.get_trajectory_index("supersaturation")].time.value=time
        self.trajectory_list[self.get_trajectory_index("growth rate")].time.value=time
        self.trajectory_list[self.get_trajectory_index("nucleation rate")].time.value=time


        C=moment_trajectory[:,10]
        u3s=moment_trajectory[:,3]
        u3n=moment_trajectory[:,8]
        for i in range(len(C)):
            T=self.trajectory_list[self.get_trajectory_index("temperature")].trajectory.value[i]
            Csat=self.solubility(T)
            S=(C[i]-Csat)/Csat
            if S<0:
                S=0
            G=self.kg.value[0]*S**self.g.value[0]
            B=self.kb.value[0]*S**self.b.value[0]*(u3s[i]+u3n[i])
            
            self.trajectory_list[self.get_trajectory_index("supersaturation")].trajectory.value.append(S)
            self.trajectory_list[self.get_trajectory_index("growth rate")].trajectory.value.append(G)
            self.trajectory_list[self.get_trajectory_index("nucleation rate")].trajectory.value.append(B)


        # update tau and f0 trajectory
        self.trajectory_list[self.get_trajectory_index("transform time")].time.value=time
        self.trajectory_list[self.get_trajectory_index("PDF for size=0")].time.value=time

        self.trajectory_list[self.get_trajectory_index("transform time")].trajectory.value.append(0)
        self.trajectory_list[self.get_trajectory_index("PDF for size=0")].trajectory.value.append(0)

        for i in range(len(time)-1):

            # update tau
            delta_t=time[i+1]-time[i]
            G_mean=(self.trajectory_list[self.get_trajectory_index("growth rate")].trajectory.value[i+1]+self.trajectory_list[self.get_trajectory_index("growth rate")].trajectory.value[i])/2
            tau_old=self.trajectory_list[self.get_trajectory_index("transform time")].trajectory.value[i]
            self.trajectory_list[self.get_trajectory_index("transform time")].trajectory.value.append(tau_old+G_mean*delta_t)
            
            # update f0
            G_ins=self.trajectory_list[self.get_trajectory_index("growth rate")].trajectory.value[i+1]
            B_ins=self.trajectory_list[self.get_trajectory_index("nucleation rate")].trajectory.value[i+1]
            if G_ins==0:
                self.trajectory_list[self.get_trajectory_index("PDF for size=0")].trajectory.value.append(0)
            else:
                self.trajectory_list[self.get_trajectory_index("PDF for size=0")].trajectory.value.append(B_ins/G_ins)

        # Update PDF
        size=[]
        PDF=[]
        
        tau_end=self.trajectory_list[self.get_trajectory_index("transform time")].trajectory.value[-1]
        tmp1=len(self.trajectory_list[self.get_trajectory_index("transform time")].trajectory.value)
        for i in range(tmp1):
            index=tmp1-1-i
            size.append(tau_end-self.trajectory_list[self.get_trajectory_index("transform time")].trajectory.value[index])
            PDF.append(self.trajectory_list[self.get_trajectory_index("PDF for size=0")].trajectory.value[index])

        seed_size=list_operation(self.seed_PDF.size.value,"+",tau_end)
        seed_PDF=self.seed_PDF.density_function.value

        size.extend(seed_size)
        PDF.extend(seed_PDF)

        self.product_PDF.size.value=size
        self.product_PDF.density_function.value=PDF


        # print the calculation time
        end_time=tm.time()
        program_time=round(end_time-start_time,4)
        string="execution time= "+str(program_time)+" s"
        print(string)

    def output(self,setting_dir,CSD_mesh_size):
        """
        Export the simulation result (PDF, CSD and trajectory)
        """

        # export the PDF
        PDF_output=[self.product_PDF.size,self.product_PDF.density_function]
        variable_output(PDF_output,setting_dir+"Result/PDF.csv")

        # export the CSD
        mesh_size=CSD_mesh_size
        size,number_fraction=DF_to_NF(self.product_PDF.size,self.product_PDF.density_function,mesh_size)
        CSD_output=[size,number_fraction]
        variable_output(CSD_output,setting_dir+"Result/CSD(number).csv")

        size,volume_fraction=NF_to_VF(size,number_fraction)
        CSD_output=[size,volume_fraction]
        variable_output(CSD_output,setting_dir+"Result/CSD(volume).csv")


        # export the trajectory
        trajectory_output=[]
        trajectory_output.append(self.trajectory_list[0].time)
        for i in range(len(self.trajectory_list)):
            trajectory_output.append(self.trajectory_list[i].trajectory)
        variable_output(trajectory_output,setting_dir+"Result/trajectory.csv")


        



    


    
    

        
    