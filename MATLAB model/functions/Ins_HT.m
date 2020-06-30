function [ Q ] = Ins_HT( T,rho,T_w,V,dV_dT,w,d_piston,gamma,heat_transfer)
%% This function calculates the heat transfer from cylinder wall to the refrigerant
if heat_transfer == 1

%%%%%%%%
%Approximate Viscosity and Conductivity
%%%%%%%%
mu=5.87509602E-07+3.79308232E-08*T;  %approximation for mu between 290 and 390K, in kg/m-s
k = -2.218332E-02 + 1.690639E-04*T + -1.560088E-07*T^2; %Curve fit for conductivity for R134 between 290K and 390K @ 500 kPa, in W/m-K
Cp = refpropm('C','T',T,'D',rho,'R134a');

D_h = d_piston/100;
A_ht = (V/D_h)*pi*D_h;
A_piston = pi*(D_h/2)^2;

Pr = (mu*Cp)/k;
u = abs(0.5*dV_dT*(w/A_piston));
Re = (rho*u*D_h)/mu;

%%%%%%%%%%%%%%%
% Adair correlation
%%%%%%%%%%%%%%%

h_c = 0.053*(k/D_h)*Pr^(0.6)*Re^(0.8);

%%%%%%%%%%%%%%
% Instant Heat Transfer
%%%%%%%%%%%%%%

Q = h_c*A_ht*(T_w-T);       %Net heat into control volume

elseif heat_transfer == 0
    
    Q=0;
    
else
    error('Incorrect input for heat_transfer, only 0 or 1')
end

end