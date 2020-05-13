function [a] = x_RK_flux_dom(x_1,x_2,rho,A_valve,k_valve,V,m_eff,C_d,A_port)
%% this fucntion calculates the derivative of valve lift in flux dominant region

a=[x_2;(rho*A_port*(V-x_2)^2-k_valve*x_1+0.5*C_d*rho*V^2*A_valve-0.5*C_d*rho*x_2^2*A_valve)/m_eff];


end

