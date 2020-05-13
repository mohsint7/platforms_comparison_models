function [ m_dot_in,m_dot_out,x_valve_suc,x_dot_valve_suc,x_valve_dis,x_dot_valve_dis ] = valve1(P_i,P_d,P,rho_i,rho,T_i,T,R,gamma,valve_dynamics,x_valve_suc_1,x_dot_valve_suc_1,x_valve_dis_1,x_dot_valve_dis_1,step,w)
%Valve flow model


%Calculate Valve Port Areas
d_discharge=0.0059;                 %discharge port diameter in meters


d_suction=d_discharge;               %suction diameter in meters

% valve flow area
A_discharge=pi*((d_discharge^2)/4);
A_suction=(pi*((d_suction^2)/4));

A_port = A_discharge;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This model is for 'digital valves'
% This means that the valves are either fully open or fully closed
% This model uses isentropic nozzel flow model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if valve_dynamics == 0
        if P_i>P    %suction valve is open
            PR=P/P_i;     "ratio of discharge pressure to the pressure inside the cylinder";
            PR_c=(1+(gamma-1)/2)^(gamma/(1-gamma));   "Critical pressure ratio";
            rho_up=P_i*1000.0/(R*T_i);
    if PR>PR_c % chocked flow
         m_dot_in=A_discharge*P_i*1000.0/(R*T_i)^0.5*(2*gamma/(gamma-1.0)*PR^(2.0/gamma)*(1-PR^((gamma-1.0)/gamma)))^0.5;
         m_dot_out=0;
    else 
        m_dot_in=A_discharge*rho_up*(gamma*R*T_i)^0.5*(1.+(gamma-1.)/2.)^((1+gamma)/(2*(1-gamma)));
        m_dot_out=0;

    end
        
        
        
        elseif P>P_d   % discharge valve is open
            PR=P_d/P;     "ratio of discharge pressure to the pressure inside the cylinder";
            PR_c=(1+(gamma-1)/2)^(gamma/(1-gamma));   "Critical pressure ratio";
            rho_up=P*1000.0/(R*T);
if PR>PR_c             % chocked flow
     m_dot_out=A_discharge*P*1000.0/(R*T)^0.5*(2*gamma/(gamma-1.0)*PR^(2.0/gamma)*(1-PR^((gamma-1.0)/gamma)))^0.5;
     m_dot_in=0;
else 
    m_dot_out=A_discharge*rho_up*(gamma*R*T)^0.5*(1.+(gamma-1.)/2.)^((1+gamma)/(2*(1-gamma)));
    m_dot_in=0;
  
end
        else   % both valves are closed
            m_dot_in=0;
            m_dot_out=0;
        end
% 
%     %If there is flowrate coming in...
%     if P_i>P      
%         Pr=P/P_i;
%         Ma=sqrt((2/(gamma-1))*((Pr^((1-gamma)/gamma))-1));
%        %Above calculate the Mach number based on the compressible flow
%        %relations
%         if Ma >1
%             Ma=1;
%         end
%         
%         V=sqrt(gamma*R*T_i*1000)*Ma;    %Calculate gas velocity
%         m_dot_in=rho_i*V*A_suction;     %Calculate mass flow rate
%         m_dot_out = 0;
% 
%     %if there is flowrate coming out...
%     elseif P>=P_d
%        Pr=P_d/P;
%        Ma=sqrt((2/(gamma-1))*((Pr^((1-gamma)/gamma))-1));
%        %Above calculate the Mach number based on the compressible flow
%        %relations
%        if Ma>1
%            Ma=1;
%        end
%        
%        V=sqrt(gamma*R*T*1000)*Ma;       %Gas velocity
%        m_dot_out=rho*V*A_discharge;     %mass flow rate
%        m_dot_in = 0;
%     else
%        m_dot_in=0;
%        m_dot_out = 0;
% 
%     end
%     
    %Valve motion is zero with these types of valves    
    x_valve_suc=0;
    x_dot_valve_suc=0;
    x_valve_dis=0;
    x_dot_valve_dis=0;
    
    
%%%%%%%%%%%%%%
% When dynamics are on the valves display some dynamics
%%%%%%%%%%%%%%    
elseif valve_dynamics == 1
    
    %%%%%%%%%%%%%
    % Additional Valve Givens and Calculations
    %%%%%%%%%%%%%
        
    d_valve = 0.0059;         %Valve diameter, meters
    A_valve = pi*(d_valve/2)^2;
    C_d = 1.17;             %Drag coefficient
    E = 1.93e11;            %Youngs Modulus, Pa
    h_valve = 0.0001532;    %Valve thickness, meters
    l_valve = 0.018;            %total length of valve, meters
    a_valve = 0.0140;            %distance from anchor to force, meters
    rho_valve = 8000;           %density of spring steel, kg/m^3
    dis_stopper = 0.0018;        %discharge stopper height, meters
    
    x_tr=0.25*(d_discharge^2/d_valve);  %Transitional lift
    
    I=(d_valve*h_valve^3)/12;  %Moment of Intertia for discharge valve,m^4
    
    k_valve=(6*E*I)/(a_valve^2*(3*l_valve-a_valve));    %Valve stiffness
    m_eff=(1/3)*rho_valve*l_valve*d_valve*h_valve;      %Effective mass of valve reeds 
    
    %If there is flowrate coming in...
        if P_i>P
            PR=P/P_i;     "ratio of discharge pressure to the pressure inside the cylinder";
            PR_c=(1+(gamma-1)/2)^(gamma/(1-gamma));   "Critical pressure ratio";
            rho_up=P_i*1000.0/(R*T_i);
    if PR>PR_c
         m_dot_in=A_discharge*P_i*1000.0/(R*T_i)^0.5*(2*gamma/(gamma-1.0)*PR^(2.0/gamma)*(1-PR^((gamma-1.0)/gamma)))^0.5;
         m_dot_out=0;
    else 
        m_dot_in=A_discharge*rho_up*(gamma*R*T_i)^0.5*(1.+(gamma-1.)/2.)^((1+gamma)/(2*(1-gamma)));
        m_dot_out=0;

    end
        
        V=m_dot_in/(rho_up*A_valve); 
        
        P_high = P_i;
        P_low = P;
        
        %Mass-flux dominant
        if x_valve_suc_1>=x_tr

            [a] = x_RK_flux_dom(x_valve_suc_1,x_dot_valve_suc_1,rho_i,A_valve,k_valve,V,m_eff,C_d,A_port);

        %Pressure dominant
        else 
            
            [a] = x_RK_pressure_dom(x_valve_suc_1,x_dot_valve_suc_1,rho_i,A_valve,k_valve,V,m_eff,C_d,P_high,P_low);


        end

        x_valve_suc = x_valve_suc_1 + a(1)*(step/w);
        x_dot_valve_suc = x_dot_valve_suc_1 + a(2)*(step/w);
        
        %Logic that keeps the valve from moving beyond the stopper
        if x_valve_suc >dis_stopper
            x_valve_suc = dis_stopper;
            x_dot_valve_suc = 0;
        end
        
        A_massflow = pi*d_valve*x_valve_suc;
         m_dot_in=A_massflow*P_i*1000.0/(R*T_i)^0.5*(2*gamma/(gamma-1.0)*PR^(2.0/gamma)*(1-PR^((gamma-1.0)/gamma)))^0.5;
         m_dot_out=0;
        x_valve_dis = 0;
        x_dot_valve_dis = 0;
        
          
    %if there is flowrate coming out...
    elseif P>P_d
            PR=P_d/P;     "ratio of discharge pressure to the pressure inside the cylinder";
            PR_c=(1+(gamma-1)/2)^(gamma/(1-gamma));   "Critical pressure ratio";
            rho_up=P*1000.0/(R*T);
if PR>PR_c
     m_dot_out=A_discharge*P*1000.0/(R*T)^0.5*(2*gamma/(gamma-1.0)*PR^(2.0/gamma)*(1-PR^((gamma-1.0)/gamma)))^0.5;
     m_dot_in=0;
else 
    m_dot_out=A_discharge*rho_up*(gamma*R*T)^0.5*(1.+(gamma-1.)/2.)^((1+gamma)/(2*(1-gamma)));
    m_dot_in=0;
  
end
       
       V=m_dot_out/(rho_up*A_valve);
       
       P_high = P;
       P_low = P_d;
    
        %Mass flux dominant 
        if x_valve_dis_1>=x_tr

            [a] = x_RK_flux_dom(x_valve_dis_1,x_dot_valve_dis_1,rho,A_valve,k_valve,V,m_eff,C_d,A_port);
            
        %Pressure dominant
        else 
            
            [a] = x_RK_pressure_dom(x_valve_dis_1,x_dot_valve_dis_1,rho,A_valve,k_valve,V,m_eff,C_d,P_high,P_low);


        end
        
        x_valve_dis = x_valve_dis_1 + a(1)*(step/w);
        x_dot_valve_dis = x_dot_valve_dis_1 + a(2)*(step/w);
        
        %Logic that keeps the valve from moving beyond the stopper
        if x_valve_dis >dis_stopper
            x_valve_dis = dis_stopper;
            x_dot_valve_dis = 0;
        end

        A_massflow = pi*d_valve*x_valve_dis;
        m_dot_out=A_massflow*P*1000.0/(R*T)^0.5*(2*gamma/(gamma-1.0)*PR^(2.0/gamma)*(1-PR^((gamma-1.0)/gamma)))^0.5;
     m_dot_in=0;
        x_valve_suc = 0;
        x_dot_valve_suc = 0;
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   %If the pressure is between the suction and discharge pressure
   %there is no flow in or out
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   else
       m_dot_in=0;
       m_dot_out = 0;
       x_valve_suc = 0;
       x_dot_valve_suc = 0;
       x_valve_dis = 0;
       x_dot_valve_dis = 0;

   end
    
    
else
    
    error('Incorrect Input for valve_dynamics flag')
    
end


end

