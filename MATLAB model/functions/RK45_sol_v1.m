function x=RK45_sol(dtheta,V_dead,V_disp,theta,rhoi,Ti,du_drho,du_dT,w,hin,P_s,P_d,T0,rho0,k,R,d,Qdot,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,error_allowed)
%% A simplified function based on RK45 integration technique which will be used for cycle integration in compressor model
g=1;
%initialization
error=1;           
error_max=1;

%error tolerances
%error_allowed=1e-8;
%dtheta_min=1e-4;
while error_max>error_allowed;
t(1)=theta;
rho(1)=rhoi; 
T(1)=Ti;

P1=refpropm('P','T',T(1),'D',rho(1),'R134a');               %kPa
[mdot_in1,mdot_out1,x_valve_suc1,x_dot_valve_suc1,x_valve_dis1,x_dot_valve_dis1 ] = valve1(P_s,P_d,P1,rho0,rho(1),T0,T(1),R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
y1=ind_prop_der(V_dead,V_disp,t(1),rho(1),T(1),du_drho,du_dT,w,hin,mdot_in1,mdot_out1,Qdot);

%% f2
t(2)=t(1)+(1/5)*dtheta; 
rho(2)=rho(1)+(1/5)*dtheta*y1(1);   T(2)=T(1)+(1/5)*dtheta*y1(2);
P2=refpropm('P','T',T(2),'D',rho(2),'R134a');               %kPa
[mdot_in2,mdot_out2,x_valve_suc2,x_dot_valve_suc2,x_valve_dis2,x_dot_valve_dis2 ] = valve1(P_s,P_d,P2,rho0,rho(2),T0,T(2),R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
y2=ind_prop_der(V_dead,V_disp,t(2),rho(2),T(2),du_drho,du_dT,w,hin,mdot_in2,mdot_out2,Qdot);

%% f3
t(3)=t(1)+(3/10)*dtheta;  rho(3)=rho(1)+dtheta*((3/40)*y1(1)+(9/40)*y2(1));   T(3)=T(1)+dtheta*((3/40)*y1(2)+(9/40)*y2(2));
P3=refpropm('P','T',T(3),'D',rho(3),'R134a');               %kPa
[mdot_in3,mdot_out3,x_valve_suc3,x_dot_valve_suc3,x_valve_dis3,x_dot_valve_dis3 ] = valve1(P_s,P_d,P3,rho0,rho(3),T0,T(3),R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
y3=ind_prop_der(V_dead,V_disp,t(3),rho(3),T(3),du_drho,du_dT,w,hin,mdot_in3,mdot_out3,Qdot);

%% f4
t(4)=t(1)+(3/5)*dtheta;  rho(4)=rho(1)+dtheta*((3/10)*y1(1)-(9/10)*y2(1)+(6/5)*y3(1));   T(4)=T(1)+dtheta*((3/10)*y1(2)-(9/10)*y2(2)+(6/5)*y3(2));
P4=refpropm('P','T',T(4),'D',rho(4),'R134a');               %kPa
[mdot_in4,mdot_out4,x_valve_suc4,x_dot_valve_suc4,x_valve_dis4,x_dot_valve_dis4 ] = valve1(P_s,P_d,P4,rho0,rho(4),T0,T(4),R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
y4=ind_prop_der(V_dead,V_disp,t(4),rho(4),T(4),du_drho,du_dT,w,hin,mdot_in4,mdot_out4,Qdot);

%% f5
t(5)=t(1)+dtheta;  rho(5)=rho(1)+dtheta*(-(11/54)*y1(1)+(5/2)*y2(1)-(70/27)*y3(1)+(35/27)*y4(1));   T(5)=T(1)+dtheta*(-(11/54)*y1(2)+(5/2)*y2(2)-(70/27)*y3(2)+(35/27)*y4(2));
P5=refpropm('P','T',T(5),'D',rho(5),'R134a');              %kPa
[mdot_in5,mdot_out5,x_valve_suc5,x_dot_valve_suc5,x_valve_dis5,x_dot_valve_dis5 ] = valve1(P_s,P_d,P5,rho0,rho(5),T0,T(5),R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
y5=ind_prop_der(V_dead,V_disp,t(5),rho(5),T(5),du_drho,du_dT,w,hin,mdot_in5,mdot_out5,Qdot);

%% f6
t(6)=t(1)+(7/8)*dtheta;  rho(6)=rho(1)+dtheta*((1631/55296)*y1(1)+(175/512)*y2(1)+(575/13824)*y3(1)+(44275/110592)*y4(1)+(253/4096)*y5(1));  
T(6)=T(1)+dtheta*((1631/55296)*y1(2)+(175/512)*y2(2)+(575/13824)*y3(2)+(44275/110592)*y4(2)+(253/4096)*y5(2));
P6=refpropm('P','T',T(6),'D',rho(6),'R134a');               %kPa
[mdot_in6,mdot_out6,x_valve_suc6,x_dot_valve_suc6,x_valve_dis6,x_dot_valve_dis6 ] = valve1(P_s,P_d,P6,rho0,rho(6),T0,T(6),R,k,valve_dynamics,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis,dtheta,w);
y6=ind_prop_der(V_dead,V_disp,t(6),rho(6),T(6),du_drho,du_dT,w,hin,mdot_in6,mdot_out6,Qdot);

%% error
e_rho=dtheta*(-(277/64512)*y1(1)+(6925/370944)*y3(1)-(6925/202752)*y4(1)-(277/14336)*y5(1)+(277/7084)*y6(1));
e_T=dtheta*(-(277/64512)*y1(2)+(6925/370944)*y3(2)-(6925/202752)*y4(2)-(277/14336)*y5(2)+(277/7084)*y6(2));
error=[e_rho,e_T];
error_max=max(abs(error));


% take smaller step and reiterate if condition is met
if error_max>error_allowed
dtheta=0.9*dtheta*(error_allowed/error_max)^0.3;
end
g=g+1;
end
theta1=theta+dtheta;
% next state of the refrigerant
rho=rho(1)+dtheta*((37/378)*y1(1)+(250/621)*y3(1)+(125/594)*y4(1)+(512/1771)*y6(1));
T=T(1)+dtheta*((37/378)*y1(2)+(250/621)*y3(2)+(125/594)*y4(2)+(512/1771)*y6(2));

% take bigger step for next step
if error_max<error_allowed
dtheta1=0.9*dtheta*(error_allowed/error_max)^0.2;
end


x=[rho,T,theta1,dtheta,x_valve_suc1,x_dot_valve_suc1,x_valve_dis1,x_dot_valve_dis1,mdot_in1,mdot_out1,dtheta1];
end