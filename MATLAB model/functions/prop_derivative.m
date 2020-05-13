function [du_dT,du_drho]= prop_derivative(T,rho)
T_1=T;
rho_1=rho;

dT=0.1;
drho=0.1;

T_2=T_1+dT;
rho_2=rho_1+drho;


u_T2=refpropm('U','T',T_2,'D',rho_1,'R134a');
u_T1=refpropm('U','T',T_1,'D',rho_1,'R134a');

u_rho1=refpropm('U','T',T_1,'D',rho_1,'R134a');
u_rho2=refpropm('U','T',T_1,'D',rho_2,'R134a');

du_dT=(u_T2-u_T1)/(T_2-T_1);

du_drho=(u_rho2-u_rho1)/(rho_2-rho_1);
end



