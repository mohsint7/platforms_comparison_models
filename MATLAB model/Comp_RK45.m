 clear all
clc
Folder = cd;
addpath('functions');
addpath('C:\Program Files (x86)\REFPROP');
format long
tic
%% number of steps in one cycle
solver_tol=1e-8;             % solver tolerance
dtheta0=1e-7;                %initial step size to start the solver solver
valve_dynamics = input('Turn on valve dynamics? 1 for on 0 for off: ');  %Zero for off, One for on
heat_transfer = input('Turn on heat transfer? 1 for on 0 for off: ');    %Zero for off, One for on

%% Compressor parameters

Vdead=8e-8;          %clearance volume [m]
V_disp=8e-6;         % Displacement volume [m]
d=0.0059;             %suction and discharge valve diameter in m%


B=2;                 %cylinder bore diameter in cm%
w=377;               % Angular speed [rad/s]
N=60*w/(2*pi);       % RPM
PR=2.5 ;              %compressor pressure ratio





%% initial conditions & suction side states
rho0=23.75;                                          %density,[kg/m3], R134a%
T0=293;                                             %eveaporation temperature or compressure inlet temperature[K]%

%% fluid properties
h_in = refpropm('H','T',T0,'D',rho0,'R134a');               %suction side enthalpy [J/kg]
u=refpropm('U','T',T0,'D',rho0,'R134a');                    %J/kg
P=refpropm('P','T',T0,'D',rho0,'R134a');                    %kP

R=0.08314*1000;                                             %specific gas constant[J/kg.k]


P_s=refpropm('P','T',T0,'D',rho0,'R134a');                      %Kpa
P_d=P_s*PR;                                                    %Kpa




%% Initializing vectors for vlave model

x_valve_suc(1)=0;
x_dot_valve_suc(1)=0;
x_valve_dis(1)=0;
x_dot_valve_dis(1)=0;


%Property derivatives for compression process equations
[du_dT,du_drho]= prop_derivative(T0,rho0);



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

%% Mass and Energy balance
% vectors initialization
 f=1 ;
 
  error=1;
  rad(1)=0;
  T_w(1) = 300;
   T(1)=T0;
rho(1)=rho0;
dtheta1(1)=dtheta0;
dtheta(1)=dtheta0;

while error>1e-4
i=1;
if f>1
    % resseting the vectors at after each iteration
    T_error=[];
    rho_error=[];
    T_error=T;
    rho_error=rho;
    T=[];
    rho=[];
    h=[];
    V=[];
    dV_dtheta=[];
    P=[];
    rad=[];
    k=[];
    V=[];
    dtheta=[];
    dtheta1=[];
    Qdot=[];
    Qdot=[];
    x_valve_suc=[];
    x_dot_valve_suc=[];
    x_valve_dis=[];
    x_dot_valve_dis=[];
    mdot_in=[];
    mdot_out=[];
    
    
    dtheta1(1)=dtheta0;
    dtheta(1)=dtheta0;
    x_valve_suc(1)=0;
    x_dot_valve_suc(1)=0;
    x_valve_dis(1)=0;
    x_dot_valve_dis(1)=0;
    rad(1)=0;
    
    % final state= initial state for convergence
    T(1)=T_error(end);
    rho(1)=rho_error(end);
end
while rad<2*pi
k(i)=refpropm('K','T',T0,'D',rho0,'R134a');
h(i)=refpropm('H','T',T(i),'D',rho(i),'R134a'); 
P(i)=refpropm('P','T',T(i),'D',rho(i),'R134a');               %kPa
[V(i),dV_dtheta(i)]=Volume(Vdead,V_disp,rad(i));
% instantaneous heat transfer from cylinder wall to refrigerant
[Qdot(i)]  = Ins_HT( T(i),rho(i),T_w(f),V(i),dV_dtheta(i),w,B,k(i),heat_transfer);
% Calling RK45 solver function for compression process equations
x23=RK45_sol_v1(dtheta1(i),Vdead,V_disp,rad(i),rho(i),T(i),du_drho,du_dT,w,h_in,P_s,P_d,T0,rho0,k(i),R,d,Qdot(i),valve_dynamics,x_valve_suc(i),x_dot_valve_suc(i),x_valve_dis(i),x_dot_valve_dis(i),solver_tol);
% seprating variables from solver results
rho(i+1)=x23(1);
T(i+1)=x23(2);
rad(i+1)=x23(3);
dtheta(i+1)=x23(4);
x_valve_suc(i+1)=x23(5);
x_dot_valve_suc(i+1)=x23(6);
x_valve_dis(i+1)=x23(7);
x_dot_valve_dis(i+1)=x23(8);
mdot_in(i)=x23(9);
mdot_out(i)=x23(10);
dtheta1(i+1)=x23(11);

i=i+1;


end

 P(i)=refpropm('P','T',T(i),'D',rho(i),'R134a');
 
 % Calculating volume for just plotting purposes
 [V(i),dV_dtheta(i)]=Volume(Vdead,V_disp,rad(i));
 [Qdot(i)]  = Ins_HT( T(i),rho(i),T_w(f),V(i),dV_dtheta(i),w,B,k(i-1),heat_transfer);
 % total heat transfer calculation for one cycle
Qdot1=(Qdot.*dtheta)/w;
Q_dot_cyl(f) = trapz(Qdot1);

% Friction Model

    mu_oil = 0.486;             %oil viscosity, Pa-sec
    delta_gap = 0.000050;       %Gap width, meters
    l_piston = 0.02;            %Length of piston, meters
    A_length = pi*(B/100)*l_piston;
    u_ave = 4.8;
    F_viscous = (mu_oil*A_length*u_ave)/delta_gap;
    W_dot_friction = (F_viscous*u_ave)/1000; 
    
    % cylinder wall calculation from heat transfer between ambient and
    % cylinder wall
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
% residual calculations at the end of each iteration
if f>1
res_T(f)=1-abs(max([T(1),T(end)]./[T_error(1),T_error(end)]));
res_rho(f)=1-abs(max([rho(1),rho(end)]./[rho_error(1),rho_error(end)]));
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

% total mdot just for plots
mdot=mdot_in-mdot_out;


% setting the length of the vectors for plots
rad(end)=[];
T(end)=[];
rho(end)=[];
V(end)=[];
P(end)=[];
Qdot(end)=[];
dtheta(end)=[];
dV_dtheta(end)=[];

% Compressor performance parameters calculations
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
eta_vol=m_dot_tot_in/(rho0*V_disp*(w/(2*pi)));        % volumetric efficiency

W_PV=trapz((((P*1000).*(dV_dtheta)).*dtheta)*(377/(2*pi))); % indicared power
%% Plots
subplot(3,3,1);
plot(rad,T,'k');
title (' temperature');
toc
subplot(3,3,2);
plot(rad,rho,'k');
title ('density');

subplot(3,3,3);
plot(rad,mdot,'k');
title ('mdot');
% 
% subplot(3,3,4);
% plot(theta_01,dV_dtheta,'k');
% title ('change of volume');
% 
% subplot(3,3,5);
% plot(theta_01,drho_dt,'k');
% title ('change in density');
% 
% subplot(3,3,6);
% plot(theta_01,dT_dtheta,'k');
% title ('change in temperature');

subplot(3,3,7);
plot(rad,V,'k');
title ('Volume');

subplot(3,3,8);
plot(rad,P,'k');
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
plot(rad,Qdot);title('Heat Transfer');

%% writting results to excell file 
% optional section. correct the directory address before running this
% portion of the code


Tab=table(rad',P',V',T',rho',mdot');
col_header={'theta','Pressure','Volume','Temperature','Density','Mass'};
output_matrix=[{' '} col_header ];
filename = 'C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\compressor_model_work\Software comparison work\results\PV_mat3.xlsx';

writetable(Tab,filename,'Sheet',1,'Range','B1');
xlswrite(filename,output_matrix);