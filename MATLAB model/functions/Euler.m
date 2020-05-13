function y=Euler(dtheta,V_dead,V_disp,theta,rho,T,du_drho,du_dT,w,hin,mdot_in,mdot_out,Qdot)

dy=ind_prop_der(V_dead,V_disp,theta,rho,T,du_drho,du_dT,w,hin,mdot_in,mdot_out,Qdot);

rho1=rho+dtheta*dy(1);
T1=T+dtheta*dy(2);
y=[rho1,T1];
end