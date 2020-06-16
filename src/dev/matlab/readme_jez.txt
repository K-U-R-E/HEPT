e_d = 0.1524            % Airframe internal diameter
i_d = 0.147             %Internal Motor Diameter
t_l = 0.7               % Thrust Chamber Length

% Environment

g = 9.81                %"Gravitational Acceleration, m/s^2"
Pa = 110000             %"Ambient Pressure, Pascal"
Ta = 300                %"Ambient Temperature, K

% Combustion chamber

n_cstar = 0.9           %"Combustion Efficiency"
n_cf = 0.885            %"Cf Efficiency"
k = 1.114               %"Ratio of Specific Heats with respect to combustion products"
Dint = 0.05             %"Initial Paraffin Grain Diameter, m"
Dchamber = i_d          %"Chamber Diameter, m"
pcomb = 900             %"Paraffin's density, kg/m^3"
Vol = t_l * pi * (i_d/2)^2     %"Free chamber volume, m^3"
M = 26.912              %"Molar mass of combustion products, kg/kmol"
To = 3347               %"Chamber temperature, K"
Dthroat = 26*10^-3     %"Throat Diameter, m"
Dexit = 63.36*10^-3    %"Exit Diameter, m"
Lgrain = 0.4            %"Paraffin Grain Lenght, m"

% Oxidiser Tank

Cd = 0.65               %"Injector's Discharge Coefficient"
IA = 0.0002             %"Injector Area, m^2"
dox = e_d               % Diameter of oxidiser tank (m)

Cd_Vent = 0.65          % Vent's Discharge Coefficient
i_ox_mass = 22.3        % Initial oxidiser mass, kg
length_tank =1.7        % Length of oxidiser tank, m

