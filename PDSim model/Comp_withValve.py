""""
(c) Ian Bell, 2012
www.BellThermal.com

2012 Purdue University Compressor Conference Short Course

A very simplified example of a reciprocating compressor without valves
or heat transfer to chamber

For illustrative purposes only!
"""

##########################################
###          Part 1. Imports           ###
##########################################

#Here we import some python packages
from __future__ import division, print_function
from math import pi, cos, sin
from time import clock
import os, sys
import pandas as pd
import xlsxwriter

# If the following line is uncommented, python will try to use a local version
# of PDSim.  This is handy for debugging purposes.  Generally you want this line 
# commented out
# PDSim should also be built using a command like python build_ext --inplace to keep all the extension modules next to the .pyx files
#sys.path.insert(0, os.path.abspath('..'))

#Here we import the things from PDSim we need
from PDSim.flow.flow import FlowPath
from PDSim.flow import flow_models
from PDSim.misc.datatypes import arraym
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.core import PDSimCore
from PDSim.plot.plots import debug_plots
from PDSim.flow.flow_models import ValveModel
#Imports from CoolProp (fluid property database)
from CoolProp import State
from CoolProp import CoolProp as CP
    
##########################################
###   Part 2. Declaration of PURecip   ###
##########################################

class PURecip(PDSimCore):
    def __init__(self):
        #Initialize the base class that PURecip is derived from
        PDSimCore.__init__(self)
        
    def V_dV(self, theta):
        """
        The simplest volume relationship possible, given by dead volume 
        and displacement directly
        """
        V = self.Vdead+self.Vdisp/2*(1-cos(theta))
        dVdtheta = self.Vdisp/2*sin(theta)
        return V, dVdtheta
    
    def Suction(self,FlowPath,**kwargs):
        if FlowPath.key_up=='A':
            # pressure in compressor higher than the inlet line
            # valve is closed - no flow
            return 0.0
        else:
            #try:
                
            if not hasattr(self,'valve_model') or self.valve_model  == '1DOF':
            
                FlowPath.A=self.suction_valve.A()
                mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            
            elif self.valve_model  == 'VibrationModeModel':

                #mdot = self.suction_valve.mdot
                A1=self.suction_valve.A_1()
                mdot1=flow_models.IsentropicNozzle(A1,FlowPath.State_up,FlowPath.State_down)
                A2=self.suction_valve.A_2()
                mdot2=flow_models.IsentropicNozzle(A2,FlowPath.State_up,FlowPath.State_down)   
                mdot = mdot1+mdot2         
            return mdot
        
    def Discharge(self,FlowPath,**kwargs):
        if FlowPath.key_down=='A':
            # pressure in compressor lower than the discharge line
            # valve is closed - no flow
            return 0.0
        else:
            #try:
                
            if not hasattr(self,'valve_model') or self.valve_model  == '1DOF':
            
                FlowPath.A=self.discharge_valve.A()
                mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            
            elif self.valve_model  == 'VibrationModeModel':

                #mdot = self.discharge_valve.mdot
                A1=self.discharge_valve.A_1()
                mdot1=flow_models.IsentropicNozzle(A1,FlowPath.State_up,FlowPath.State_down)
                A2=self.discharge_valve.A_2()
                mdot2=flow_models.IsentropicNozzle(A2,FlowPath.State_up,FlowPath.State_down)   
                mdot = mdot1+mdot2                
            return mdot
            
    def TubeCode(self, Tube):
        """ 
        A thin wrapper of the isothermal wall tube from flow_models.py 
        """
        Tube.Q = flow_models.IsothermalWallTube(Tube.mdot, 
                                                Tube.State1, 
                                                Tube.State2,
                                                Tube.fixed, 
                                                Tube.L, 
                                                Tube.ID,
                                                T_wall = self.Tlumps[0])
    
    def heat_transfer_callback(self, theta):
        """
        A callback used by PDSimCore.derivs to calculate the heat transfer
        to the gas in the working chamber.
        
        We return an arraym instance the same length as the number of CV in existence
        
        More code (a proper heat transfer model) could be included here, like 
        in PDSim.recip.core
        """
        return arraym([0.0]*self.CVs.N)
        
    def mechanical_losses(self):
        """
        The mechanical losses in kW
        """
        return self.Wdot_parasitic
    
    def ambient_heat_transfer(self):
        """
        The ambient heat transfer for the compressor in kW
        
        Returns a positive value if heat is added to the compressor from the 
        ambient
        """
        return self.h_shell*self.A_shell*(self.Tamb-self.Tlumps[0]) #[kW]
        
    def lump_energy_balance_callback(self):
        """
        A callback used in PDSimCore.solve to do the energy balance on the lump
        
        Note: we neglect heat transfer to the gas in the working chamber
        """
        #Mechanical losses are added to the lump
        self.Wdot_mechanical = self.mechanical_losses() #[kW]
        #Heat transfer between the shell and the ambient
        self.Qamb = self.ambient_heat_transfer() #[kW]
        return self.Wdot_mechanical + self.Qamb
 
##########################################
###    Part 3. Execution of PURecip    ###
##########################################
        
def Compressor(**kwargs):

    recip=PURecip()                     #Instantiate the model
    recip.Vdead = 8e-8                #Dead volume [m3]
    recip.Vdisp = 8e-6                  #Displacement/rev [m3]
    recip.omega = 377                  #Frequency, rad/sec (60Hz)

    recip.d_discharge=0.0059;           #discharge port diameter [m]
    recip.d_suction=recip.d_discharge;  #suction diameter [m]
    recip.A_discharge=pi*recip.d_discharge**2/4
    recip.A_suction=pi*recip.d_suction**2/4
    
    recip.EulerN= 7000
     
    recip.RK45_eps=1e-7
        
    #These are parameters needed for the ambient heat transfer model
    recip.h_shell = 0.010               #[kW/m2/K]
    recip.A_shell = pi*10*2*(0.0254**2) #[m2]
    recip.Tamb = 293                   #[K] 
    recip.plot_names='recip_plot_buttons'
    
    #Parameters for the mechanical losses model (simplified)
    recip.Wdot_parasitic = 0         #Parasitic losses [kW]
    
    Ref='R134a'
    inletState=State.State(Ref,dict(T=293.15,D=23.75))
    outletState=State.State(Ref,{'T':400,'P':inletState.p*2.5})
    mdot_guess = inletState.rho*recip.Vdisp*recip.omega/(2*pi)
    
    #First add the control volumes.
    recip.add_CV( ControlVolume(key='A',
                               initialState=outletState.copy(),
                               VdVFcn=recip.V_dV,
                               becomes='A') )
    
    #Add the inlet tube
    recip.add_tube( Tube(key1='inlet.1',key2='inlet.2',L=0.03,ID=0.01,
                         mdot=mdot_guess, State1=inletState.copy(),
                         fixed=1,TubeFcn=recip.TubeCode) )
    
    #Add the outlet tube
    recip.add_tube( Tube(key1='outlet.1',key2='outlet.2',L=0.03,ID=0.01,
                         mdot=mdot_guess, State2=outletState.copy(),
                         fixed=2,TubeFcn=recip.TubeCode) )
    
    #Add the flow paths that link flow nodes together
    recip.add_flow(FlowPath(key1='inlet.2',key2='A',MdotFcn=recip.Suction))
    recip.add_flow(FlowPath(key1='outlet.1',key2='A',MdotFcn=recip.Discharge))
    
    
    
    E = 1.93e11             #Youngs Modulus, [Pa]
    h_valve = 0.0001532     #Valve thickness, [m]
    l_valve = 0.018         #Total length of valve, [m]
    a_valve = 0.0140        #Distance from anchor to force, [m]
    rho_valve = 8000        #Density of spring steel, [kg/m^3] 
    C_D = 1.17              #Drag coefficient [-]
    d_valve = 0.0059         #Valve Diameter [m]
    x_stopper = 0.0018      #Stopper location [m]

    I=(d_valve*h_valve**3)/12  #Moment of Intertia for valve,[m^4]
    k_valve=(6*E*I)/(a_valve**2*(3*l_valve-a_valve))    #Valve stiffness
    m_eff=(1/3)*rho_valve*l_valve*d_valve*h_valve      #Effective mass of valve reeds
    x_tr_suction = 0.25*(recip.d_suction**2/d_valve)
    x_tr_discharge = 0.25*(recip.d_discharge**2/d_valve)

    #The suction valve parameters
    recip.suction_valve=ValveModel(
          d_valve=d_valve,
          d_port=recip.d_suction,
          C_D=C_D,
          rho_valve=rho_valve,
          x_stopper=x_stopper,
          m_eff = m_eff,
          k_valve = k_valve,
          x_tr = x_tr_suction,
          key_up='inlet.2',
          key_down='A'
          )
    recip.add_valve(recip.suction_valve)
    
    #The discharge valve parameters
    recip.discharge_valve=ValveModel(
          d_valve=d_valve,
          d_port=recip.d_discharge,
          C_D=C_D,
          rho_valve=rho_valve,
          x_stopper=x_stopper,
          m_eff = m_eff,
          k_valve = k_valve,
          x_tr = x_tr_discharge,
          key_up='A',
          key_down='outlet.1'
          )
    recip.add_valve(recip.discharge_valve)
    
    recip.connect_callbacks(endcycle_callback=recip.endcycle_callback, # Provided by PDSimCore
                            heat_transfer_callback=recip.heat_transfer_callback,
                            lumps_energy_balance_callback = recip.lump_energy_balance_callback
                            )
    
    t1=clock()
    recip.precond_solve(key_inlet='inlet.1',
                        key_outlet='outlet.2',
                        solver_method = 'RK45',
                        OneCycle = False,
                        UseNR = True,
                        )
    print('time taken', clock()-t1)
    debug_plots(recip)
    del recip.FlowStorage
    # from PDSim.plot.plots import debug_plots
    # debug_plots(recip)
    return recip
def run_all():
    Compressor(solver_method = 'Euler')
    Compressor(solver_method = 'Heun')
    Compressor(solver_method = 'RK45')

if __name__=='__main__':    
    #If this file is run directly, this code will be run
    recip=Compressor(solver_method = 'RK45')
    # run_all()
import matplotlib.pyplot as plt
plt.plot(recip.t,recip.xValves[0,:], lw = 1.5)
plt.show()
plt.figure()
plt.plot(recip.t,recip.xValves[2,:], lw = 1.5)
plt.show()



## exporting data to excel files
# #
# df1= pd.DataFrame(recip.V.T,recip.p.T)
# df2= pd.DataFrame(recip.t.T,recip.rho.T)
# df3= pd.DataFrame(recip.T.T,recip.m.T)
# df4=pd.DataFrame(recip.FlowsProcessed.summed_mdot['outlet.1'],recip.FlowsProcessed.summed_mdot['inlet.2'])
# df5=pd.DataFrame(recip.h.T)
# # # #
# # #writer = pd.ExcelWriter('D:\Phd\compressor_model_work\compressor_model_work\Python_codes_rec\test1.xlsx', engine='xlsxwriter')
# df1.to_excel(excel_writer = "C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\compressor_model_work\Software comparison work\model/test1.xlsx", sheet_name='Sheet1',startcol=3)
# df2.to_excel(excel_writer = "C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\compressor_model_work\Software comparison work\model/test2.xlsx", sheet_name='Sheet1',startcol=3)
# df3.to_excel(excel_writer = "C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\compressor_model_work\Software comparison work\model/test3.xlsx", sheet_name='Sheet1',startcol=3)
# df4.to_excel(excel_writer = "C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\compressor_model_work\Software comparison work\model/test4.xlsx", sheet_name='Sheet1',startcol=3)
# df5.to_excel(excel_writer = "C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\compressor_model_work\Software comparison work\model/test5.xlsx", sheet_name='Sheet1',startcol=3)