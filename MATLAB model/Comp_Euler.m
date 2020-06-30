%Simplified recriprocating compressor model with a Euler Integrator. Part
%of an electronic code annex for the work presented in Tanveer and Bradshaw
%(2020) "Quantitative and Qualitative Evaluation of Various
%Positive-Displacement Compressor Modeling Platforms" presented in Int. J.
%of Ref.
%
%This model requires an installation of REFPROP and presumes it is
%installed in the default directory. 
%
%Please reach out to mohsin.tanveer@okstate.edu or
%craig.bradshaw@okstate.edu for questions about this model.


%% Declarations and pre-processing
clear all
clc
Folder = cd;
addpath('functions');
addpath('C:\Program Files (x86)\REFPROP');

%% Model Inputs - Numerics/Flags
n=15000;            %Number of steps
tol_inner=1e-4;    %Convergence tolerancr T, rho, T_w etc
valve_dynamics = input('Turn on valve dynamics? 1 for on 0 for off: ');  %Zero for off, One for on
heat_transfer = input('Turn on heat transfer? 1 for on 0 for off: ');    %Zero for off, One for on

%% Model Inputs - Comp parameters
Vdead=8e-8;                                  %Clearance volume of the compressor
V_disp=8e-6;                                 %Displacement volume of the compressor
d=0.0059;                                    %Valve diameter in m%   
theta_01=linspace(0,360,n);                  %Crank angle%
rad=linspace(0,2*pi,n);                      %Crank angle in radian
N=3600;                                      %Compressor RPM%
B=2;                                         %Cylinder bore diameter in cm%
w=2*pi*N/60;                                 %Angular speed
PR=2.5;                                      %Compressor pressure ratio

%% Model Inputs - Fluid properties
rho0=23.75;                                  %Density,[kg/m3], R134a%
T0=293;                                      %Eveaporation temperature or compressure inlet temperature[K]%
R=81.49;                                     %Specific gas constant[J/kg.k]

%Calculation of input property data
h_in = refpropm('H','T',T0,'D',rho0,'R134a');               %Specific enthalpy at inlet [J/kg]
u=refpropm('U','T',T0,'D',rho0,'R134a');                    %Instantanious secific internal energy [J/kg]
P=refpropm('P','T',T0,'D',rho0,'R134a');                    %instantanious pressure [kPa]
P_s=refpropm('P','T',T0,'D',rho0,'R134a');                  %suction pressure [Kpa]
P_d=P_s*PR;                                                 %Dischrge side pressure [Kpa]

%% isentropic compression for isentropic eficiency calculation
%Inlet Entropy
s_1=refpropm('S','T',T0,'D',rho0,'R134a');
%Secant Method to Find Exit Temp for Isentropic Compression
dT = 1;          %Initalize change in Temperature
T_isen(1) = 325; %Guess Values
T_isen(2) = 315;
g=2;
%Secant Loop
while abs(dT)>1e-6
    T_isen(g+1)=T_isen(g)-(refpropm('S','T',T_isen(g),'P',P_d,'R134a')-s_1)*(T_isen(g)-T_isen(g-1))/(refpropm('S','T',T_isen(g),'P',P_d,'R134a')-refpropm('S','T',T_isen(g-1),'P',P_d,'R134a'));
    
    dT=T_isen(g+1)-T_isen(g);
    g=g+1;
    if g>500
        error('stuck in interations')
    end
    %keyboard
end

%Last Temperature in Iteration is Temperature at Isentropic Input
T_s_2 = T_isen(g);
%Using Temperature and Discharge Pressure, enthalpy is Calcualted
h_2_s = refpropm('H','T',T_s_2,'P',P_d,'R134a');
h_2s=refpropm('H','P',P_d,'S',s_1,'R134a');
%% Control Volume calculation
dtheta=rad(2); %fixed step size

%Property derivative for compression process equation
[du_dT,du_drho]= prop_derivative(T0,rho0);
     
%% Mass and Energy balance

%initializations
f=1;
error=1;
T_w(1) = 300;
T(1)=T0;
rho(1)=rho0;
x_valve_suc(1)=0;
x_dot_valve_suc(1)=0;
x_valve_dis(1)=0;
x_dot_valve_dis(1)=0;

%Residual tolerance loop
while error>tol_inner

    if f>1
        T_error=T;
        rho_error=rho;
        T(1)=T(n+1);
        rho(1)=rho(n+1);
    end

    %Inner loop, iterating on crank angle
    for i=1:n
        k(i)=refpropm('K','T',T(i),'D',rho(i),'R134a');
        P(i)=refpropm('P','T',T(i),'D',rho(i),'R134a');               %kPa

        % mass flow and valve model
        [mdot_in(i),mdot_out(i),x_valve_suc(i+1),x_dot_valve_suc(i+1),x_valve_dis(i+1),x_dot_valve_dis(i+1) ] = valve1(P_s,P_d,P(i),rho0,rho(i),T0,T(i),R,k(1),valve_dynamics,x_valve_suc(i),x_dot_valve_suc(i),x_valve_dis(i),x_dot_valve_dis(i),dtheta,w);
        [V(i),dV_dtheta(i)]=Volume(Vdead,V_disp,rad(i));
        % heat transfer from cylinder wall to the refrigerant
        [Qdot(i)]  = Ins_HT( T(i),rho(i),T_w(f),V(i),dV_dtheta(i),w,B,k(i),heat_transfer);       

        % Euler method function for solution
        x23=Euler(dtheta,Vdead,V_disp,rad(i),rho(i),T(i),du_drho,du_dT,w,h_in,mdot_in(i),mdot_out(i),Qdot(i));
        rho(i+1)=x23(1);
        T(i+1)=x23(2);
        i=i+1;

    end
    
    %Total heat transfer for one cycle
    Q_dot_cyl(f) = trapz(Qdot)*(dtheta/w);

    % Friction Model
    mu_oil = 0.486;             %Oil viscosity, Pa-sec
    delta_gap = 0.000050;       %Gap width, meters
    l_piston = 0.02;            %Length of piston, meters
    A_length = pi*(B/100)*l_piston;
    u_ave = 4.8;
    F_viscous = (mu_oil*A_length*u_ave)/delta_gap;
    W_dot_friction = (F_viscous*u_ave)/1000 

    %Cylinder wall temperature calculation using heat transfer to ambient  
    if heat_transfer == 1
        [Q_dot_out(f),T_w(f+1)] = outer_HT(T_w(f));
        res_HT(f) = abs(Q_dot_out(f) - Q_dot_cyl(f) - W_dot_friction);

        if f>2
            T_w(f+1) = T_w(f) - res_HT(f)*((T_w(f) - T_w(f-1))/(res_HT(f) - res_HT(f-1)));
        end
    else
        T_w(f+1) = T_w(f);
        res_HT(f) = 0;
    end

    if f>1
        %Residuals calculations
        res_T(f)=1-abs(max(T./T_error));
        res_rho(f)=1-abs(max(rho./rho_error));
        disp(res_T(f))
        disp(res_rho(f))
        error=[abs(res_T(f)),abs(res_rho(f))];
        error=max(error);
    end
    
    f=f+1

    if f>30
        error=0;
    end
end

%% Post processing --  Compressor performance parameters calculation
% total mdot just for plots
mdot=mdot_in-mdot_out;

dtime=dtheta/w;
time(1)=0;
for s=1:1:length(dtheta)
    time(s+1)=time(s)+dtime(s);
    s=s+1;
end
time(end)=[];
m_dot_tot_out =(N/60)* trapz(mdot_out.*dtime)         % average discjarge mass flow rate
m_dot_tot_in = (N/60)*trapz(mdot_in.*dtime)           % average suction mass flow rate
Wdot=m_dot_tot_out*(h_2_s-h_in)                       % isentropic power 
eta_vol=m_dot_tot_in/(rho0*V_disp*(w/(2*pi)))        % volumetric efficiency

W_PV=trapz((((P*1000).*(dV_dtheta)).*dtheta)*(377/(2*pi))) % indicared power
%% Plots
T(n+1)=[];
rho(n+1)=[];
 
subplot(3,3,1);
plot(theta_01,T,'k');
title (' temperature');

% subplot(3,3,2);
% plot(theta_01,rho,'k');
% title ('density');
% 
subplot(3,3,3);
plot(theta_01,mdot,'k');
title ('mdot');
% 
subplot(3,3,4);
plot(theta_01,dV_dtheta,'k');
title ('change of volume');
% 
% subplot(3,3,5);
% plot(theta_01,drho_dt,'k');
% title ('change in density');
% 
% subplot(3,3,6);
% plot(theta_01,dT_dtheta,'k');
% title ('change in temperature');

subplot(3,3,7);
plot(theta_01,V,'k');
title ('Volume');

subplot(3,3,8);
plot(theta_01,P,'k');
title ('Pressure');

% subplot(3,3,9);
% plot(theta_01,h,'k');
% title ('enthalpy');

figure
plot(V,P);
xlabel('Volume');
ylabel('Pressure');

figure
plot(res_T,'-r*','DisplayName','T-res');
hold on
plot(res_rho,'-.m+','DisplayName','rho-res');
xlabel('iteration');
ylabel('residual');
legend

figure
plot(theta_01,Qdot);title('Heat Transfer');

%% Exporting the results to Excel
Tab=table(rad',P',V',T',rho',mdot');
col_header={'theta','Pressure','Volume','Temperature','Density','Mass'};
output_matrix=[{' '} col_header ];
%% *******Update below to a generic filename
filename = 'PV_mat1.xlsx';

writetable(Tab,filename,'Sheet',1,'Range','B1');
xlswrite(filename,output_matrix);