"""
The code performs a simulation of a nitrous oxide - paraffin hybrid
rocket motor. Some inputs from the user are needed, as seen below. For
results, just copy it and run it at Jupyter notebook, preferably.

Author: Guilherme Tavares

Note: Additions made by Vinay Williams at Kingston University, 2020.

ASSUMPTIONS

- Isentropic flow along the nozzle

- Temperature inside the combustion chamber does not change along the
  burn

- Combustion chamber is adiabatic

- Oxidiser tank empties adiabatically (and temperature does not change in the vapour-only phase)

- Available volume of combustion chamber does not change

- Combustion products form a mixture which behaves like an ideal gas

- Ideal burning (no erosion along the burning port)

- No O/F ratio change along the burn (strong assumption, but with N20 the performance variation is not so dependant on O/F)

KNOW BEFORE RUNNING

- Some inputs should be provided by the user on the section below.

- Scipy and CoolProp needed.

- No prior computation about the ullage for the oxidiser tank is carried out by this version of the code. Carefully compute how much         oxidiser will be pumped inside the tank. Also, carefully compute the tank's diameter and length.

- The vent of the oxidiser tank is considered to be on the top of the tank.

- For plotting different variables, proceed to the bottom of the code and follow the examples.

- Step can be reduced for further precision, but must be accompanied by a reduction on "while" conditions.

- New versions of the code will be gradually uploaded to tackle the above-mentioned assumptions, and to improve user experience.

REFERENCES

For general internal ballistic's equations:
SUTTON, G.P.; BIBLARZ, Oscar. Rocket Propulsion Elements. 8th ed. John Wiley & Sons, Inc., Hoboken, New Jersey, 2010.

For self-pressurizing N20 Model:
WHITMORE, S.A.; CHANDLER, S.N. Engineering Model for Self-Pressurizing Saturated-N20-Propellant Feed Systems. Journal of Propulsion and Power, Vol. 26, No. 4, July-August 2010.

For N20-paraffin's regression law:
LESTRADE, J.-Y.; ANTHOINE, J.; LAVERGNE, G. Liquefying fuel regression rate modeling in hybrid propulsion. The French Aerospace Lab, Mauzac, France; The French Aerospace Lab, Toulouse, France, 2010.

For N20 properties:
CoolProp (available documentation in http://www.coolprop.org/index.html)

"""

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.stats import linregress
from scipy.integrate import odeint
try:
    import CoolProp.CoolProp
    from CoolProp.CoolProp import PropsSI
except ImportError:
    print('Unable to load CoolProp. CoolProp files will not be imported.')

"FLAGS"

t_s = True

"INPUTS"
e_d = 0.1524 # Airframe internal diameter
i_d = 0.147 #Internal Motor Diameter
t_l = 0.7 # Thrust Chamber Length
# Environment
g = 9.81    #"Gravitational Acceleration, m/s^2"
Pa = 110000 #"Ambient Pressure, Pascal"
Ta = 300    #"Ambient Temperature, K

# Combustion chamber
n_cstar = 0.9     #"Combustion Efficiency"
n_cf = 0.885           #"Cf Efficiency"
k = 1.114           #"Ratio of Specific Heats with respect to combustion products"
Dint = 0.05          #"Initial Paraffin Grain Diameter, m"
Dchamber = i_d     #"Chamber Diameter, m"
pcomb = 900         #"Paraffin's density, kg/m^3"
Vol = t_l * np.pi * (i_d/2) ** 2     #"Free chamber volume, m^3"
M = 26.912            #"Molar mass of combustion products, kg/kmol"
To = 3347             #"Chamber temperature, K"
Dthroat = 26*10**-3   #"Throat Diameter, m"
Dexit = 63.36*10**-3 #"Exit Diameter, m"
Lgrain = 0.4      #"Paraffin Grain Lenght, m"

# Oxidiser Tank
Cd = 0.65     #"Injector's Discharge Coefficient"
IA = 0.0002  #"Injector Area, m^2"
dox = e_d   # Diameter of oxidiser tank (m)

Cd_Vent = 0.65       # Vent's Discharge Coefficient
i_ox_mass = 22.3    # Initial oxidiser mass, kg
length_tank =1.7# Length of oxidiser tank, m

"COMPUTED PARAMETERS"

R = 8314/M; #"Universal Gas Constant divided by Effective Molecular Mass"

VA = np.pi*(0.001**2)/4 # Vent Area, m2
gamma_vapour =  PropsSI('C','T',Ta,'Q',0,'NITROUSOXIDE')/PropsSI('CVMASS','T',Ta,'Q',0,'NITROUSOXIDE') # Specific heat ratio of vapour nitrous oxide
Vtank = np.pi*dox**2*length_tank/4; # Volume of oxidiser tank, m3

X = (PropsSI('D','T',Ta,'Q',1,'NITROUSOXIDE')*PropsSI('D','T',Ta,'Q',0,'NITROUSOXIDE')*Vtank-PropsSI('D','T',Ta,'Q',1,'NITROUSOXIDE')*i_ox_mass)/(i_ox_mass*(PropsSI('D','T',Ta,'Q',0,'NITROUSOXIDE')
-PropsSI('D','T',Ta,'Q',1,'NITROUSOXIDE'))) # Initial quality of nitrous oxide inside the tank
if X < 0:
    print("Warning: X is less than zero! Increase tank volume or decrease initial mass")
    exit()

"INITIALIZE MONITOR SOLUTIONS"

mass_nitrous = [i_ox_mass] # Available nitrous mass at each time instant, kg
mass_liquid_nitrous = [(1-X)*i_ox_mass] # Available liquid nitrous mass at each time instant, kg
mass_vapour_nitrous = [X*i_ox_mass] # Available vapour nitrous mass at each time instant, kg
temperature_tank = [Ta] # Tank temperature at each time instant, K
out_liquid_nitrous = [0] # Mass of liquid nitrous oxide that escapes through the injector, kg
out_vapour_nitrous = [0] # Mass of vapour nitrous oxide that escapes through the injector, kg
consumed_liquid = [0] # Total consumed liquid nitrous at each time instant, kg
consumed_vapour = [0] # Total consumed vapour nitrous at each time instant, kg
oxy = [0]
chamber_pressure = [Pa] # Combustion chamber pressure, Pa
grain_radius = [Dint/2] # Paraffin grain radius, m
tank_pressure = [PropsSI('P','T',Ta,'Q',X,'NITROUSOXIDE')] # Oxidiser tank pressure, Pa
total_entropy = [i_ox_mass*PropsSI('S','T',Ta,'Q',X,'NITROUSOXIDE')] # Total entropy with respect to oxidiser tank, J/K
specific_entropy = [PropsSI('S','T',Ta,'Q',X,'NITROUSOXIDE')] # Specific entropy with respect to oxidiser tank, J/kgK
overall_density = [i_ox_mass/Vtank] # Density of nitrous oxide inside the tank (considering both phases), kg/m3
quality = [X] # Quality of the satured nitrous oxide
thrust = [0] # Thrust generated by the motor, N
cf = [0] # Thrust coefficient

########################################
# TWO-PHASE FLOW (INSIDE OXIDISER TANK)
########################################

def two_phase(x,t,k,M,R,To,Lgrain,Cd,IA,pcomb,n_cstar,Dthroat,Vol,T):
        "x[0] is the chamber pressure Po, Pascal"
        "x[1] is the radius of Paraffin Grain, m"
        "x[2] is the used liquid oxidiser mass, kg"
        "x[3] is the used vapour oxidiser mass, kg"
        "T is the temperature inside the tank, in K"
        pho_liquid_nitrous = PropsSI('D','T',T,'Q',0,'NITROUSOXIDE')
        vapour_pressure_nitrous = PropsSI('P','T',T,'Q',1,'NITROUSOXIDE')
        gamma_vapour =  PropsSI('C','T',T,'Q',1,'NITROUSOXIDE')/PropsSI('CVMASS','T',T,'Q',1,'NITROUSOXIDE')
        oxidiser_flow_rate = Cd*IA*(2*pho_liquid_nitrous*(vapour_pressure_nitrous-x[0]))**0.5
        oxy.append(oxidiser_flow_rate)
        pgas = x[0]*M/(8314*To);
        dx0 = (R*To/Vol)*(oxidiser_flow_rate+Lgrain*2*np.pi*x[1]*(pcomb-pgas)*1.28*10**(-5)*((oxidiser_flow_rate/(np.pi*(x[1]**2)))**0.94)-x[0]/n_cstar*(np.pi*(Dthroat**2)/4)*((k/(R*To))*(2/(k+1))**((k+1)/(k-1)))**0.5)
        dx1 = 1.28*10**(-5)*(oxidiser_flow_rate/(np.pi*(x[1]**2)))**0.94
        dx2 = oxidiser_flow_rate
        dx3 = VA*Cd_Vent*(gamma_vapour*PropsSI('P','T',T,'Q',1,'NITROUSOXIDE')*PropsSI('D','T',T,'Q',1,'NITROUSOXIDE')*(2/(gamma_vapour+1))**( (gamma_vapour+1)/(gamma_vapour-1) ))**0.5
        dx = [dx0, dx1, dx2, dx3]
        return dx


step = 0.01
t = []

x0 = [chamber_pressure[0],grain_radius[0],0,0]

ysol = []

i = 0
T = Ta
print("Two Phase Calculations....\n")
while mass_liquid_nitrous[i] > 0.01:
    ysol.append(x0)
    ts = [step*i,step*(i+1)]
    y = odeint(two_phase,x0,ts,args=(k,M,R,To,Lgrain,Cd,IA,pcomb,n_cstar,Dthroat,Vol,T))
    x0 = y[1,:].tolist()
    t.append(i*step)
    i = i + 1
    out_liquid_nitrous.append(x0[2]-consumed_liquid[i-1])
    out_vapour_nitrous.append(x0[3]-consumed_vapour[i-1])
    total_entropy.append(total_entropy[i-1]-out_liquid_nitrous[i]*PropsSI('S','T',T,'Q',0,'NITROUSOXIDE')-out_vapour_nitrous[i]*PropsSI('S','T',T,'Q',1,'NITROUSOXIDE'))
    mass_nitrous.append(mass_nitrous[i-1]-out_liquid_nitrous[i]-out_vapour_nitrous[i])
    specific_entropy.append(total_entropy[i]/mass_nitrous[i])
    overall_density.append(mass_nitrous[i]/Vtank)
    T = PropsSI('T','S',specific_entropy[i],'D',overall_density[i],'NITROUSOXIDE')
    temperature_tank.append(T)
    quality.append(PropsSI('Q','S',specific_entropy[i],'D',overall_density[i],'NITROUSOXIDE'))
    mass_liquid_nitrous.append((1-quality[i])*mass_nitrous[i])
    mass_vapour_nitrous.append(quality[i]*mass_nitrous[i])
    consumed_liquid.append(x0[2])
    consumed_vapour.append(x0[3])
    chamber_pressure.append(x0[0])
    grain_radius.append(x0[1])
    tank_pressure.append(PropsSI('P','T',T,"Q",quality[i],'NITROUSOXIDE'))
    cf.append((((2*k**2)/(k-1))*(2/(k+1))**((k+1)/(k-1))*(1-(Pa/chamber_pressure[i])**((k-1)/k)))**0.5)
    thrust.append(n_cf*cf[i]*(np.pi*(Dthroat**2)/4)*chamber_pressure[i])

#########################################
# VAPOUR-ONLY PHASE (INSIDE OXIDISER TANK)
#########################################

def vapour_phase(x,t,k,M,R,To,Lgrain,Cd,IA,pcomb,n_cstar,Dthroat,Vol,P):
        "x[0] is the chamber pressure Po, Pascal"
        "x[1] is the radius of Paraffin Grain, m"
        "x[2] is the consumed vapour oxidiser mass, kg"
        "P is the pressure inside the tank, in Pascal"
        pho_vapour_nitrous = PropsSI('D','P',P,'Q',1,'NITROUSOXIDE')
        gamma_vapour =  PropsSI('C','P',P,'Q',1,'NITROUSOXIDE')/PropsSI('CVMASS','P',P,'Q',1,'NITROUSOXIDE')
        oxidiser_flow_rate = Cd*IA*(gamma_vapour*PropsSI('P','T',T,'Q',1,'NITROUSOXIDE')*PropsSI('D','T',T,'Q',1,'NITROUSOXIDE')*(2/(gamma_vapour+1))**( (gamma_vapour+1)/(gamma_vapour-1) ))**0.5
        oxy.append(oxidiser_flow_rate)
        pgas = x[0]*M/(8314*To);
        dx0 = (R*To/Vol)*(oxidiser_flow_rate+Lgrain*2*np.pi*x[1]*(pcomb-pgas)*1.28*10**(-5)*((oxidiser_flow_rate/(np.pi*(x[1]**2)))**0.94)-x[0]/n_cstar*(np.pi*(Dthroat**2)/4)*((k/(R*To))*(2/(k+1))**((k+1)/(k-1)))**0.5)
        dx1 = 1.28*10**(-5)*(oxidiser_flow_rate/(np.pi*(x[1]**2)))**0.94
        dx2 = oxidiser_flow_rate + VA*Cd_Vent*(gamma_vapour*PropsSI('P','T',T,'Q',1,'NITROUSOXIDE')*PropsSI('D','T',T,'Q',1,'NITROUSOXIDE')*(2/(gamma_vapour+1))**( (gamma_vapour+1)/(gamma_vapour-1) ))**0.5
        dx = [dx0, dx1, dx2]
        return dx

x0 = [chamber_pressure[i],grain_radius[i],consumed_vapour[i]]
ysol = []

P = PropsSI('P','T',temperature_tank[i],"Q",quality[i],'NITROUSOXIDE')
print("Vapour Phase Calculations....\n")
while mass_vapour_nitrous[i] > 0.01 and P>Pa:
    ysol.append(x0)
    ts = [step*i,step*(i+1)]
    y = odeint(vapour_phase,x0,ts,args=(k,M,R,To,Lgrain,Cd,IA,pcomb,n_cstar,Dthroat,Vol,P))
    x0 = y[1,:].tolist()
    t.append(i*step)
    i = i + 1
    out_liquid_nitrous.append(0)
    out_vapour_nitrous.append(x0[2]-consumed_vapour[i-1])
    total_entropy.append(total_entropy[i-1]-out_vapour_nitrous[i]*PropsSI('S','P',P,'Q',1,'NITROUSOXIDE'))
    mass_nitrous.append(mass_nitrous[i-1]-out_vapour_nitrous[i])
    specific_entropy.append(total_entropy[i]/mass_nitrous[i])
    overall_density.append(mass_nitrous[i]/Vtank)
    #T = PropsSI('T','S',specific_entropy[i],'D',overall_density[i],'NITROUSOXIDE')
    temperature_tank.append(temperature_tank[i-1])
    P = PropsSI('P','T',temperature_tank[i-1],'D',overall_density[i],'NITROUSOXIDE')
    quality.append(1)
    mass_liquid_nitrous.append(mass_liquid_nitrous[i-1])
    mass_vapour_nitrous.append(mass_nitrous[i])
    consumed_liquid.append(consumed_liquid[i-1])
    consumed_vapour.append(x0[2])
    chamber_pressure.append(x0[0])
    grain_radius.append(x0[1])
    tank_pressure.append(P)
    cf.append((((2*k**2)/(k-1))*(2/(k+1))**((k+1)/(k-1))*(1-(Pa/chamber_pressure[i])**((k-1)/k)))**0.5)
    thrust.append(n_cf*cf[i]*(np.pi*(Dthroat**2)/4)*chamber_pressure[i])

################
# TAIL-OFF PHASE
################

def tail_off(x,t,k,M,R,To,n_cstar,Dthroat,Vol):
        "x is the chamber pressure Po, Pascal"
        pgas = x*M/(8314*To);
        dx0 = (R*To/Vol)*(-x/n_cstar*(np.pi*(Dthroat**2)/4)*((k/(R*To))*(2/(k+1))**((k+1)/(k-1)))**0.5)
        return dx0
print("Tail off Calculations....\n")
ysol = []
x0 = chamber_pressure[i]
while chamber_pressure[i] > 1.1*Pa:
    ysol.append(x0)
    ts = [step*i,step*(i+1)]
    y = odeint(tail_off,x0,ts,args=(k,M,R,To,n_cstar,Dthroat,Vol))
    t.append(i*step)
    x0 = y[1]
    i = i + 1
    chamber_pressure.append(x0)
    total_entropy.append(total_entropy[i-1])
    specific_entropy.append(specific_entropy[i-1])
    overall_density.append(overall_density[i-1])
    temperature_tank.append(temperature_tank[i-1])
    quality.append(1)
    mass_nitrous.append(mass_nitrous[i-1])
    mass_liquid_nitrous.append(mass_liquid_nitrous[i-1])
    mass_vapour_nitrous.append(mass_vapour_nitrous[i-1])
    consumed_liquid.append(consumed_liquid[i-1])
    consumed_vapour.append(consumed_vapour[i-1])
    grain_radius.append(grain_radius[i-1])
    tank_pressure.append(tank_pressure[i-1])
    if Pa < chamber_pressure[i]:
        cf.append((((2*k**2)/(k-1))*(2/(k+1))**((k+1)/(k-1))*(1-(Pa/chamber_pressure[i])**((k-1)/k)))**0.5)
        thrust.append(n_cf*cf[i]*(np.pi*(Dthroat**2)/4)*chamber_pressure[i])
    else:
        cf.append(0)
        thrust.append(0)

######
# PLOT
######

print("Plotting Figures....\n")

fig = plt.figure()
fig.suptitle('Thrust Curve', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Thrust (N)', fontsize=10)
plt.plot(t,thrust[0:len(thrust)-1])
plt.savefig("./images/thrustcurve.png")
fig2 = plt.figure()
fig2.suptitle('Pressure Curve', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Pressure (Pa)', fontsize=10)
plt.plot(t,chamber_pressure[0:len(chamber_pressure)-1])
plt.savefig("./images/pressurecurve.png")
fig3 = plt.figure()
fig3.suptitle('Tank Pressure ', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Pressure (Pa)', fontsize=10)
plt.plot(t,tank_pressure[0:len(tank_pressure)-1])
plt.savefig("./images/tankpressure.png")
fig4 = plt.figure()
fig4.suptitle('Tank Temperature', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Temperature (K)', fontsize=10)
plt.plot(t,temperature_tank[0:len(temperature_tank)-1])
plt.savefig("./images/tanktemperature.png")
fig5 = plt.figure()
fig5.suptitle('Grain Radius ', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Diameter (Meters)', fontsize=10)
plt.plot(t,grain_radius[0:len(grain_radius)-1])
plt.savefig("./images/grainradius.png")
fig6 = plt.figure()
fig6.suptitle('Quality', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Phase Ratio', fontsize=10)
plt.plot(t,quality[0:len(quality)-1])
plt.savefig("./images/quality.png")
fig7 = plt.figure()
fig7.suptitle('Mass N2O', fontsize=12)
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Mass N2O (kg)', fontsize=10)
plt.plot(t,mass_nitrous[0:len(quality)-1])
plt.savefig("./images/n2omass.png")

import pandas as pd
print(len(t))
print(len(mass_nitrous))
print(len(thrust))
print(len(chamber_pressure))
print(len(grain_radius))
print(len(temperature_tank))
print(len(quality))

# todo@vinay: augment arrays to be the same length, that is time array is one unit smaller than the other array
# todo@vinay: because the dataframe won't work right if they aren't the same length
data = {"Time":t, "Mass":mass_nitrous, "Thrust":thrust, "Pressure":chamber_pressure, "Grain Radius":grain_radius, "Tank Temperature":temperature_tank, "Tank Pressure":tank_pressure, "Quality":quality}
df = pd.DataFrame(data, columns = ["Time", "Thrust", "Pressure", "Grain Radius", "Tank Temperature", "Tank Pressure", "Quality"])
df.to_csv("data1.csv", index = False)


slope, intercept, r_value, p_value, std_err = linregress(t,mass_nitrous[0:len(quality)-1])
print("N20 Mass Flow Rate: "+str(slope))
print("Final Grain Radius: "+str(round(grain_radius[len(mass_nitrous)-1],5)))
print(" ")

totalforce = np.sum(thrust)[0]
time = t[len(t)-1]
impulse = (totalforce/len(thrust))* time
print("Impulse: "+ str(round(impulse,5)))
