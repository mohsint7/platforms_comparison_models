function mdot=discharge_valve1(P,P_s,T,rho,k,R,d,x)
PR=P/P_s;     "ratio of discharge pressure to the pressure inside the cylinder";
PR_c=(1+(k-1)/2)^(k/(1-k));   "Critical pressure ratio";
rho_up=P_s*1000.0/(R*T);
A=(pi*d^2)/4;
if P<P_s
if PR>PR_c
     mdot=A*P*1000.0/(R*T)^0.5*(2*k/(k-1.0)*PR^(2.0/k)*(1-PR^((k-1.0)/k)))^0.5;
else 
    mdot=A*rho_up*(k*R*T)^0.5*(1.+(k-1.)/2.)^((1+k)/(2*(1-k)));
  
end
else
    mdot=0;
end