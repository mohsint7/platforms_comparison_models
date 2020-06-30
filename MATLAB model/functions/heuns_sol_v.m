function y2=heuns_sol(dtheta,V_dead,V_disp,theta,rho,T,du_drho,du_dT,w,hin,P_s,P_d,T0,rho0,k,R,d,Qdot,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis)

%% A simplified function based on Heun's method which will be used for cylce integration in compressor model

P1=refpropm('P','T',T,'D',rho,'R134a');               %kPa

[ mdot_in1,mdot_out1,x_valve_suc1,x_dot_valve_suc1,x_valve_dis1,x_dot_valve_dis1 ] = valve1(P_s,P_d,P1,rho0,rho,T0,T,R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
 

dy1=ind_prop_der(V_dead,V_disp,theta,rho,T,du_drho,du_dT,w,hin,mdot_in1,mdot_out1,Qdot);
rho1=rho+dtheta*dy1(1);
T1=T+dtheta*dy1(2);


% y2
P2=refpropm('P','T',T1,'D',rho1,'R134a');               %kPa
[ mdot_in2,mdot_out2,x_valve_suc2,x_dot_valve_suc2,x_valve_dis2,x_dot_valve_dis2 ] = valve1(P_s,P_d,P2,rho0,rho,T0,T,R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
 

dy2=ind_prop_der(V_dead,V_disp,theta,rho1,T1,du_drho,du_dT,w,hin,mdot_in2,mdot_out2,Qdot);

% mdot_in=(mdot_in1+mdot_in2)/2;
% mdot_out=(mdot_out1+mdot_out2)/2;
% x_valve_suc=(x_valve_suc1+x_valve_suc2)/2;
% x_dot_valve_suc=(x_dot_valve_suc1+x_dot_valve_suc2)/2;
% x_valve_dis=(x_valve_dis1+x_valve_dis2)/2;
% x_dot_valve_dis=(x_dot_valve_dis1+x_dot_valve_dis2)/2;
rho2=rho+dtheta/2*(dy1(1)+dy2(1));
T2=T+dtheta/2*(dy1(2)+dy2(2));


y2=[rho2,T2,x_valve_suc1,x_dot_valve_suc1,x_valve_dis1,x_dot_valve_dis1, mdot_in1,mdot_out1];

%% y_f

end