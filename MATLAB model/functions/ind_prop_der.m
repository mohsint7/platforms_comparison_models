function x=ind_prop_der(V_dead,V_disp,theta,rho,T,du_drho,du_dT,w,hin,mdot_in,mdot_out,Qdot)
%% this function returns the derivative of T and rho at a step which can be used to calculte the properties at next time step
h = refpropm('H','T',T,'D',rho,'R134a');            %J/kg
u =refpropm('U','T',T,'D',rho,'R134a');              %J/kg


[V,dV_dtheta]=Volume(V_dead,V_disp,theta);

drho_dt=(1/V)*(-rho*dV_dtheta+(1/w)*(mdot_in-mdot_out));
dT_dtheta=((-rho*h*dV_dtheta)-((u*V+rho*V*du_drho)*drho_dt)+(1/w)*(Qdot+mdot_in*hin-mdot_out*h))/(rho*V*du_dT);
x=[drho_dt,dT_dtheta];
end