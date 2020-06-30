function [Q,T_w_new] = outer_HT(T_w)
%% This function calculates the cylinder wall temperature
k_alum=0.0262;                  %Thermal conductivity of R134a, W/m-K
T_amb=22+273.15;                %Ambient Temperature

rho_air=1.2;                    %Density of air, kg/m^3
Vel_air=1;                      %Velocity of airflow across compressor, m/s
D_shell=2*0.0254;               %Outer Diameter of Piston Cylinder, m
L_shell=10*0.0254;              %Length of Compressor, m
A_shell=pi*D_shell*L_shell;     %Surface area of compressor shell, m^2

mu_air=2*10^-5;                 %Dynamic Viscosity of air, kg/m-s

nu_air=mu_air/rho_air;          %Kinematic Viscosity of air, m^2/s
alpha_air=2.2160*10^-5;         %Thermal diffusivity of air, m^2/s
Pr=nu_air/alpha_air;            %Prandtl Number of air flow
Re_air=(rho_air*Vel_air*D_shell)/mu_air;    %Reynolds number of air
C_2=0.683;                      %Constants from Table 7.2 in Incropera DeWitt
m=0.466;

Nu_D=C_2*Re_air^m*Pr^(1/3);

h_shell=(Nu_D*k_alum)/D_shell;  %Heat transfer coefficient for outside heat transfer

R_shell=1/(h_shell*A_shell);    %Overall thermal resistance from shell to air

Q = h_shell*A_shell*(T_amb - T_w);
Q= Q;

T_w_new = (Q/(h_shell*A_shell))+T_amb;