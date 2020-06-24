within ;
package Compressor_model
  package functions
    " mass flow and heat transfer functions to be used in compressor model"
  function suction_valve   "This function calculates mass flow thorugh suction valve and is based on isentropic flow through nozzel. 
  This model assummes onyl either valve is fully open or close. This is a pressure dominant model. It will only opens when inside 
 pressure is less than suction side pressure."

    /************************************* Imports *********************************************/
    import SI = Modelica.SIunits;

    /************************************* Inputs **********************************************/

    input SI.Pressure P                            "Instantaneous pressure inside the cylinder";
    input SI.Pressure P_s                          "Suction side pressure";
    input SI.Temperature T                         "Instantaneous temperature inside the cylinder";
    input SI.Density  rho                          "Instantaneous density inside the cylinder";
    input SI.RatioOfSpecificHeatCapacities gamma   "Heat capacity ratio";
    input Real R                                   "Gas constant";
    input SI.Diameter d                            "Valve port diameter";

    /******************************* Intermediate variables ************************************/

    SI.Area A_v                                    "Valve port area";
    Real M_s                                       "Suction flow mach number";
    Real PR(min=0)                                 "Pressure ratio";
    SI.Velocity C_e                                "Sound velocity at specific fluid conditions inside the cylinder";

    /************************************** Outputs *****************************************/
    output SI.MassFlowRate Mdot_s                  "Suction flow rate";


  algorithm
    // If suction pressure is higher than cylinder pressure
     if P<P_s then
    PR:=P_s/P;
    A_v:=Modelica.Constants.pi*d^2/4;
    M_s:=((PR^((gamma-1)/gamma)-1)*(2/(gamma-1)))^0.5;
    // Choked Flow condition
    if M_s>1 then
      M_s:=1;
    end if;
    C_e:=sqrt(gamma*R*T);
    Mdot_s:=0.58*rho*M_s*C_e*A_v;
     else
       // Closed valve
       Mdot_s:=0;
     end if;




  end suction_valve;

  function suction_valve1   "This function calculates mass flow thorugh suction valve and is based on isentropic flow through nozzel. 
  This model assummes onyl either valve is fully open or close. This is a pressure dominant model. It will only opens when inside 
 pressure is less than suction side pressure."

    /************************************* Imports *********************************************/
    import SI = Modelica.SIunits;

    /************************************* Inputs **********************************************/

    input SI.Pressure P                            "Instantaneous pressure inside the cylinder";
    input SI.Pressure P_s                          "Suction side pressure";
    input SI.Temperature T                         "Instantaneous temperature inside the cylinder";
    input SI.Density  rho                          "Instantaneous density inside the cylinder";
    input SI.RatioOfSpecificHeatCapacities gamma   "Heat capacity ratio";
    input Real R                                   "Gas constant";
    input SI.Diameter d                            "Valve port diameter";

    /******************************* Intermediate variables ************************************/

    SI.Area A_v                                    "Valve port area";
    Real PR(min=0)                                 "Pressure ratio";
    SI.Velocity C_e                                "Sound velocity at specific fluid conditions inside the cylinder";
    Real PR_c                                     "Critical Pressure ratio";
    SI.Density rho_up                                    "Upstream density";

    /************************************** Outputs *****************************************/
    output SI.MassFlowRate Mdot_s                  "Suction flow rate";


  algorithm
    // If suction pressure is higher than cylinder pressure
     if P_s>P then
    PR:=P/P_s;
    A_v:=Modelica.Constants.pi*d^2/4;
    PR_c:=(1 + (gamma - 1)/2)^(gamma/(1 - gamma));
    rho_up:=P_s/(R*T);
    // Choked Flow condition
    if PR>PR_c then
      Mdot_s:=A_v*P/(R*T)^0.5*(2*gamma/(gamma - 1.0)*PR^(2.0/
          gamma)*(1 - PR^((gamma - 1.0)/gamma)))^0.5;
     else
       Mdot_s:=A_v*rho_up*(gamma*R*T)^0.5*(1. + (gamma - 1.)/2.)^((1 +
          gamma)/(2*(1 - gamma)));
    end if;
     else
       Mdot_s:=0;
     end if;


  end suction_valve1;

  function suction_velocity   "This function calculates mass flow thorugh suction valve and is based on isentropic flow through nozzel. 
  This model assummes onyl either valve is fully open or close. This is a pressure dominant model. It will only opens when inside 
 pressure is less than suction side pressure."

    /************************************* Imports *********************************************/
    import SI = Modelica.SIunits;

    /************************************* Inputs **********************************************/

    input SI.Pressure P                            "Instantaneous pressure inside the cylinder";
    input SI.Pressure P_s                          "Suction side pressure";
    input SI.Temperature T                         "Instantaneous temperature inside the cylinder";
    input SI.Density  rho                          "Instantaneous density inside the cylinder";
    input SI.RatioOfSpecificHeatCapacities gamma   "Heat capacity ratio";
    input Real R                                   "Gas constant";

    /******************************* Intermediate variables ************************************/

    SI.Area A_v                                    "Valve port area";
    Real M_s                                       "Suction flow mach number";
    Real PR(min=0)                                 "Pressure ratio";
    SI.Velocity C_e                                "Sound velocity at specific fluid conditions inside the cylinder";

    /************************************** Outputs *****************************************/
    output SI.Velocity V_s                  "Suction flow rate";


  algorithm
    // If suction pressure is higher than cylinder pressure
     if P<P_s then
    PR:=P_s/P;
    M_s:=((PR^((gamma-1)/gamma)-1)*(2/(gamma-1)))^0.5;
    // Choked Flow condition
    if M_s>1 then
      M_s:=1;
    end if;
    C_e:=sqrt(gamma*R*T);
    V_s:=0.58*M_s*C_e;
     else
       // Closed valve
       V_s:=0;
     end if;




  end suction_velocity;

  function discharge_valve   "This function calculates mass flow thorugh Discahrge valve and is based on isentropic flow through nozzel. 
  This model assummes onyl either valve is fully open or close. This is a pressure dominant model. It will only opens when inside 
  pressure is higher than Discharge side pressure."

    /************************************* Imports *********************************************/
    import SI = Modelica.SIunits;


    /************************************* Inputs *********************************************/
    input SI.Pressure P                                               "Pressure inside the cylinder";
    input SI.Pressure P_d                                             "Discharge side pressure";
    input SI.Temperature T                                            "Instantaneous temperature inside the cylinder";
    input SI.Density  rho                                             "Instantaneous density inside the cylinder";
    input SI.RatioOfSpecificHeatCapacities gamma                      "Heat capacity ratio";
    input Real R                                                      "Gas constant";
    input SI.Diameter d                                               "Valve port diameter";

    /***************************** Intermediate Variables ************************************/

    SI.Area A_v                                                       "Valve area";
    Real M_d                                                          "Discharge flow mach number";
    Real PR                                                           "Pressure ratio";
    SI.Velocity C_d                                                   "Sound velocity at specific fluid conditions inside the cylinder";

    /************************************* Outputs *******************************************/
    output SI.MassFlowRate Mdot_d                                     "Discharge mass flow rate";


  algorithm
    // If Discahrge pressure is less than cylinder pressure
     if P>P_d or P==P_d then
    PR:=P/P_d;
    A_v:=Modelica.Constants.pi*d^2/4;
    M_d:=((PR^((gamma-1)/gamma)-1)*(2/(gamma-1)))^0.5;
    // Choked flow condition
    if M_d>1 then
      M_d:=1;
    end if;
    C_d:=sqrt(gamma*R*T);
    Mdot_d:=0.58*rho*M_d*C_d*A_v;
    // Vlave is close
     else
       Mdot_d:=0;
     end if;




  end discharge_valve;

  function discharge_valve1   "This function calculates mass flow thorugh Discahrge valve and is based on isentropic flow through nozzel. 
  This model assummes onyl either valve is fully open or close. This is a pressure dominant model. It will only opens when inside 
  pressure is higher than Discharge side pressure."

    /************************************* Imports *********************************************/
    import SI = Modelica.SIunits;


    /************************************* Inputs *********************************************/
    input SI.Pressure P                                               "Pressure inside the cylinder";
    input SI.Pressure P_d                                             "Discharge side pressure";
    input SI.Temperature T                                            "Instantaneous temperature inside the cylinder";
    input SI.Density  rho                                          "upstream instantaneous density";
    input SI.RatioOfSpecificHeatCapacities gamma                      "Heat capacity ratio";
    input Real R                                                      "Gas constant";
    input SI.Diameter d                                               "Valve port diameter";

    /***************************** Intermediate Variables ************************************/

    SI.Area A_v                                                       "Valve area";
    Real PR_c                                                          "critical pressure ratio";
    Real PR                                                           "Pressure ratio";
    SI.Velocity C_d                                                   "Sound velocity at specific fluid conditions inside the cylinder";
    SI.Density rho_up;

    /************************************* Outputs *******************************************/
    output SI.MassFlowRate Mdot_d                                     "Discharge mass flow rate";


  algorithm
    // If Discahrge pressure is less than cylinder pressure
     if P>P_d or P==P_d then
    PR:=P_d/P;
    A_v:=Modelica.Constants.pi*d^2/4;
    PR_c:=(1 + (gamma - 1)/2)^(gamma/(1 - gamma));
    rho_up:=P/(R*T);
    // Choked flow condition
    if PR>PR_c then
      Mdot_d:=A_v*P/(R*T)^0.5*(2*gamma/(gamma - 1.0)*PR^(2.0/gamma)*(1 -
          PR^((gamma - 1.0)/gamma)))^0.5;
    else
    Mdot_d:=A_v*rho_up*(gamma*R*T)^0.5*(1. + (gamma - 1.)/2.)^((1 +
          gamma)/(2*(1 - gamma)));
    end if;
     else
       Mdot_d:=0;
     end if;




  end discharge_valve1;

  function Inside_HT  "This function calculates the heat transfer inside the compressor control volume.
  It assumes the velocity of fluid equal to the velocity of the piston. Heat transfer coefficient is 
  being calculated using Adeir et al. For more details please refer to the report.  "

    /************************************* Imports *********************************************/
    import SI = Modelica.SIunits;

    /************************************* Inputs *********************************************/

    input SI.Temperature T             "Instantaneous temperature inside the compressor";
    input SI.Density rho               "Instantaneous Density inside the compressor";
    input SI.Temperature T_w           "Cylinder wall temperature";
    input SI.Volume V                  "Instantaneous volume inside the compressor";
    input SI.VolumeFlowRate dV_dt      "Instantaneous Volume derivartive inside the compressor";
    input SI.Length d_p                "piston diameter";
    input SI.DynamicViscosity mu       "Dynamic viscosity of refrigerant";
    input Real Pr                      "Prandtle number";
    input SI.ThermalConductivity k     "Thermal conductivity of the fluid";


    /****************************** Intermediate Variables ***********************************/

    SI.Area A_ht                       "Heat transfer area";
    SI.Area A_p                        "piston area";
    SI.CoefficientOfHeatTransfer h_c   "Convection heat transfer coefficient";
    Real Re                            "Reynold's number";
    SI.Velocity Vel                    "Velocity of the fluid";

    /************************************* Outputs *********************************************/

    output SI.HeatFlowRate Q           "Heat transfer";




  algorithm

    A_ht:=(V/d_p);
    A_p:=Modelica.Constants.pi*(d_p/2)^2;
    Vel:=abs(dV_dt/A_p);
    Re:=rho*d_p*Vel/mu;

    h_c:=0.053*(k/d_p)*Pr^0.6*Re^0.8;
    Q:=h_c*A_ht*(T_w - T);


  end Inside_HT;

  function tic  "Function to record the internal time just before the experiment starts in [s]"
    Real ms_tic;
    Real sec_tic;
    Real min_tic;
    Real hour_tic;
    Real day_tic;
    Real mon_tic;
    Real year_tic;

  output Modelica.SIunits.Time tic
  "Record the internal time at the execution of this (tic) function";

  algorithm

  (ms_tic,sec_tic,min_tic,hour_tic,day_tic,mon_tic,year_tic) :=Modelica.Utilities.System.getTime();

  tic := (ms_tic*0.001) + (sec_tic) + (min_tic*60) + (hour_tic*3600);

  Modelica.Utilities.Streams.print("tic =" + String(tic) + "[s]");


  end tic;

  function toc  "Function to record the internal time just after the experiment starts in [s]"
    Real ms_toc;
    Real sec_toc;
    Real min_toc;
    Real hour_toc;
    Real day_toc;
    Real mon_toc;
    Real year_toc;

  output Modelica.SIunits.Time toc
  "Record the internal time at the execution of this (toc) function";

  algorithm

  (ms_toc,sec_toc,min_toc,hour_toc,day_toc,mon_toc,year_toc) :=Modelica.Utilities.System.getTime();

  toc := (ms_toc*0.001) + (sec_toc) + (min_toc*60) + (hour_toc*3600);

  Modelica.Utilities.Streams.print("toc =" + String(toc) + "[s]");


  end toc;

    function outside_HT "This function calculate the heat transfer from cylinder wall to the environment. 
  It need some input for air velocity and only consider convective heat transfer to the ambient.
  It do not calcualte the readiation heat transfer. Also, it do not consider shell heat transfer for 
  hermetic compressors. For more details, please refer to the report."

      /************************************* Imports *****************************************/
      import SI = Modelica.SIunits;

      /************************************* Inputs ******************************************/
      input SI.Temperature T_amb                "Ambient temperature";
      input SI.DynamicViscosity mu_air          "Dynamic viscosity of the air";
      input SI.ThermalConductivity k_solid      "Thermal condictivity of cylinder wall";
      input SI.Density rho_air                  "Ambient Air density";
      input SI.Velocity V_Air                   "Ambient Air Velocity";
      input SI.Diameter D_shell                 "outer diameter of cylinder shel";
      input SI.Length L_shell                   "Length of the cylinder shell";
      input Real Pr_Air                         "Prandtle number of the air";

      /**************************** Intermediate Variables **********************************/

      SI.Area A_shell                           "surface area of the shell";
      Real Re                                   "Reynolds number";
      Real Nu                                   "Nustle number";
      SI.HeatFlowRate Q;
      SI.CoefficientOfHeatTransfer h            "Convection heat transfer coefficient of ambient air";
      SI.Temperature T_w                        "Wall Temperature";

      /************************************ Outputs ****************************************/
      output SI.Temperature T_new;

    algorithm

    A_shell:=Modelica.Constants.pi*D_shell*L_shell;
    Re:=(rho_air*V_Air*D_shell)/mu_air;
    Nu:=0.683*Re^0.466*Pr_Air^(1/3);
    h:=(Nu*k_solid)/D_shell;

    Q:=h*A_shell*(T_amb - T_w);

    T_new:=(Q/(h*A_shell)) + T_amb;




    end outside_HT;
  end functions;

package complete_model
  " Compressor model with heat transfer and without heat transfer. Rec_Compressor_v is incomplete valve model due to modelcia dynamics"
  model Rec_Compressor_heat
    "This is the coplete model for piston cylinder compressor"

    /****************************Imports*************************************/
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;

    /*********************Compressor Specifications******************************/
    parameter SI.Length B=(0.02)                              "cylinder bore diamter in m";
    parameter SI.Volume V_dead=8e-8                   "clearance Volume in m3";
    parameter SI.Volume V_disp=8e-6;
    parameter SI.Diameter D_shell=0.05;
    parameter SI.Length L_shell=0.254;
    parameter SI.Diameter d=0.0059                            "Valve diameters";
    parameter SI.AngularVelocity w =   from_rpm(3600);
    parameter Real PR=2.5                                       "Compression ratio";

    /******************Instantaneous Fluid Properties***************************/
    SI.Density rho                             "Initial values of density in Kg/m3";
    SI.Temperature T                           "Initial value of temperature in in K";
    SI.Pressure P(displayUnit="kPa")           "Pressure at specified temperature and density conditions";
    SI.SpecificEnthalpy h                      "Enthalpy at specified temperature and density conditions";
    SI.SpecificEnergy u                        "Specific Internal enerygy";
    SI.DynamicViscosity mu=11.668e-6           "Dynamic Viscosity";
    parameter Real Pr=0.8                      "Prandtle number";

    /******************Suction Fluid Properties***************************/

    SI.Pressure P_s                           "Suction Pressure";
    parameter SI.Temperature T_s=293          "Suction Tempersture";
    parameter SI.Density rho_s=23.75           "Suction Density";
    SI.SpecificEnthalpy h_s                   "suction flow specific Enthalpy";

    /******************Discharge Fluid Properties***************************/
    SI.Pressure P_d;

    /******************ThermoPhysical Properties***************************/
    parameter SI.RatioOfSpecificHeatCapacities gamma=1.26       "Specific heat capacity ratio";
    parameter Real R=81.49                                       "Gas Constant";
    Real drhodT=-6.1634e3;
    Real dudT=5.4693e3;

    /****************** Outside Ambient Conditions ***************************/
    parameter SI.DynamicViscosity mu_air=16.82e-6;
    parameter SI.ThermalConductivity k_solid=160   "Shell thermal conductivity";
    parameter SI.Density rho_air=1.2               "Air density";
    parameter SI.Velocity V_air=1                  "Wind velocity";
    parameter Real Pr_air=0.7                      "Outside air Prandtle number";

    /****************** Mass Flow Rates ***************************/
    SI.MassFlowRate mdot_in;
    SI.MassFlowRate mdot_out;
    SI.MassFlowRate mdot;

    /****************** Intermediate Variables ***************************/
    SI.Angle theta_01(displayUnit="Deg")                           "Crank angle in rad";
    SI.Volume V(displayUnit="cm3")                                 "Volume at each time step inside the cylinder";
    SI.VolumeFlowRate dVdt(displayUnit="cm3/s")                   "Volume derivative w.r.t time";
    SI.HeatFlowRate Qdot;
    SI.Temperature Tamb=295;
    SI.Temperature T_w;

    /****************** VLE Fluid import ***************************/

  TILMedia.VLEFluid_dT vleFluid(
    computeTransportProperties=true,
                                T=T,d=rho,
    final vleFluidType=sim.vleFluidType1)
    annotation (Placement(transformation(extent={{-12,46},{8,66}})));
  TILMedia.VLEFluid_dT vleFluid1(
    computeTransportProperties=true,
                                 T=T_s,d=rho_s,
    final vleFluidType=sim.vleFluidType1)
    annotation (Placement(transformation(extent={{-14,-2},{6,18}})));

    inner TIL.SystemInformationManager sim(redeclare
        TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType1)
      annotation (Placement(transformation(extent={{54,66},{74,86}})));

  initial equation
    T=T_s;
    rho=rho_s;

  equation

    // Fluid properties of refrigerant mixture inside the cmpressor
    P=vleFluid.p;
    h=vleFluid.h;

    u=h-P/rho;

    //suction side fluid properties

    P_s=vleFluid1.p;
    h_s=vleFluid1.h;

    // Discharge side Pressure
    P_d=P_s*PR;

    /*************************** Volume Calculation ***************************/
    der(theta_01)= w;
    V = V_dead+V_disp/2*(1-cos(theta_01));
    dVdt=der(V);

    /********************************** Mass Flow *******************************/

    mdot_in=functions.suction_valve1(P,P_s,T_s,rho_s,gamma,R,d);
    mdot_out=functions.discharge_valve1(P,P_d,T,rho,gamma,R,d);
    mdot=mdot_in-mdot_out;

    /***************************** Heat Transfer ***************************/

    Qdot=functions.Inside_HT(T,rho,T_w,V,dVdt,B,mu,Pr,13.5);
    T_w=functions.outside_HT(Tamb,mu_air,k_solid,rho_air,V_air,D_shell,L_shell,Pr_air);

    /************************* Mass Energy Balance **********************/

    der(rho)=(1/V)*(-rho*dVdt+(mdot_in-mdot_out));
    der(T)=((-rho*h*dVdt)-((u*V+rho*V*drhodT)*der(rho))+(Qdot+mdot_in*h_s-mdot_out*h))/(rho*V*dudT);

  annotation (experiment(
          StopTime=0.02,
          __Dymola_NumberOfIntervals=1000,
          Tolerance=1e-07,
          __Dymola_Algorithm="Cerk45"));
  end Rec_Compressor_heat;

model Rec_Compressor_noheat
  "This is the coplete model for piston cylinder compressor"

  /****************************Imports*************************************/
  import SI = Modelica.SIunits;
  import Modelica.SIunits.Conversions.*;

  /*********************Compressor Specifications******************************/
  parameter SI.Length B=(0.02)                              "Cylinder bore diameter";
  parameter SI.Volume V_dead=8e-8                   "Clearance volume";
  parameter SI.Volume V_disp=8e-6                   "Displacement volume";
  parameter SI.Diameter D_shell=0.05                "Piston diameter";
  parameter SI.Length L_shell=0.254                 "Length of the shell";
  parameter SI.Diameter d=0.0059                     "Valve diameter";
  parameter SI.AngularVelocity w =   from_rpm(3600)   "Angular velocity";
  parameter Real PR=2.5                                       "Compression ratio";

  /******************Instantaneous Fluid Properties***************************/
  SI.Density rho                             "Initial values of density in Kg/m3";
  SI.Temperature T                           "Initial value of temperature in in K";
  SI.Pressure P(displayUnit="kPa")           "Pressure at specified temperature and density conditions";
  SI.SpecificEnthalpy h                      "Enthalpy at specified temperature and density conditions";
  SI.SpecificEnergy u                        "Specific Internal enerygy";
  SI.DynamicViscosity mu=11.668e-6           "Dynamic Viscosity";
  parameter Real Pr=0.8                      "Prandtl number";

  /******************Suction Fluid Properties***************************/

  SI.Pressure P_s                           "Suction Pressure";
  parameter SI.Temperature T_s=293          "Suction temperature";
  parameter SI.Density rho_s=23.75          "Suction density";
  SI.SpecificEnthalpy h_s                   "suction flow specific Enthalpy";

  /******************Discharge Fluid Properties***************************/
  SI.Pressure P_d;

  /******************ThermoPhysical Properties***************************/
  parameter SI.RatioOfSpecificHeatCapacities gamma=1.26       "Specific heat capacity ratio";
  parameter Real R=81.49                                       "Gas constant";
  Real drhodT=-6.1634e3;
  Real dudT=5.4693e3;

  /****************** Outside Ambient Conditions ***************************/
  parameter SI.DynamicViscosity mu_air=16.82e-6;
  parameter SI.ThermalConductivity k_solid=160   "Shell thermal conductivity";
  parameter SI.Density rho_air=1.2               "Air density";
  parameter SI.Velocity V_air=1                  "Wind velocity";
  parameter Real Pr_air=0.7                      "Outside air Prandtl number";

  /****************** Mass Flow Rates ***************************/
  SI.MassFlowRate mdot_in;
  SI.MassFlowRate mdot_out;
  SI.MassFlowRate mdot;

  /****************** Intermediate Variables ***************************/
  SI.Angle theta_01(displayUnit="Deg")                           "Crank angle in rad";
  SI.Volume V(displayUnit="cm3")                                 "Volume at each time step inside the cylinder";
  SI.VolumeFlowRate dVdt(displayUnit="cm3/s")                   "Volume derivative w.r.t time";
  SI.HeatFlowRate Qdot;
  SI.Temperature Tamb=295;
  SI.Temperature T_w;
  Real Wdot_PV;

  Modelica.SIunits.Time T_initial;
  Modelica.SIunits.Time T_final;

  /****************** VLE Fluid import ***************************/

TILMedia.VLEFluid_dT vleFluid(
  computeTransportProperties=true,
    interpolateTransportProperties=false,
                              T=T,d=rho,
  final vleFluidType=sim.vleFluidType1)
  annotation (Placement(transformation(extent={{-12,46},{8,66}})));
TILMedia.VLEFluid_dT vleFluid1(
  computeTransportProperties=true,
    interpolateTransportProperties=false,
                               T=T_s,d=rho_s,
  final vleFluidType=sim.vleFluidType1)
  annotation (Placement(transformation(extent={{-14,-2},{6,18}})));

  inner TIL.SystemInformationManager sim(redeclare
      TILMedia.VLEFluidTypes.Refprop_R134a vleFluidType1)
    annotation (Placement(transformation(extent={{54,64},{74,84}})));

initial equation
  T=T_s;
  rho=rho_s;

algorithm
  T_initial:=functions.tic();

equation

  // Fluid properties of refrigerant mixture inside the cmpressor
  P=vleFluid.p;
  h=vleFluid.h;

  u=h-P/rho;

  //suction side fluid properties

  P_s=vleFluid1.p;
  h_s=vleFluid1.h;

  // Discharge side Pressure
  P_d=P_s*PR;

  /*************************** Volume Calculation ***************************/
  der(theta_01)= w;
  V = V_dead+V_disp/2*(1-cos(theta_01));
  dVdt=der(V);

  /********************************** Mass Flow *******************************/

  mdot_in=functions.suction_valve1(P,P_s,T_s,rho_s,gamma,R,d);
  mdot_out=functions.discharge_valve1(P,P_d,T,rho,gamma,R,d);
  mdot=mdot_in-mdot_out;

  /***************************** Heat Transfer ***************************/

  //Qdot=functions.Inside_HT(T,rho,T_w,V,dVdt,B,mu,Pr,13.5);
  Qdot=0;
  T_w=functions.outside_HT(Tamb,mu_air,k_solid,rho_air,V_air,D_shell,L_shell,Pr_air);

  /************************* Mass Energy Balance **********************/

  der(rho)=(1/V)*(-rho*dVdt+(mdot_in-mdot_out));
  der(T)=((-rho*h*dVdt)-((u*V+rho*V*drhodT)*der(rho))+(Qdot+mdot_in*h_s-mdot_out*h))/(rho*V*dudT);

  Wdot_PV=P*dVdt;

algorithm

   T_final:=functions.toc();

annotation (experiment(
      StopTime=0.0173,
      __Dymola_NumberOfIntervals=1000,
      __Dymola_Algorithm="Cerk45"));
end Rec_Compressor_noheat;

model Rec_Compressor_v
  "This is the coplete model for piston cylinder compressor"

  /****************************Imports*************************************/
  import SI = Modelica.SIunits;
  import Modelica.SIunits.Conversions.*;

  /*********************Compressor Specifications******************************/
  parameter SI.Length R_01=(0.0115)                         "Crank radius in m";
  parameter SI.Length R_12=(0.1055)                         "connecting rod length";
  parameter SI.Length B=(0.02)                              "cylinder bore diamter in m";
  parameter SI.Volume V_dead=(0.00000038)                   "clearance Volume in m3";
  parameter SI.Diameter D_shell=0.05;
  parameter SI.Length L_shell=0.254;
  parameter SI.Diameter d=0.0059                            "Valve diameters";
  parameter SI.AngularVelocity w =   from_rpm(3000);
  parameter Real PR=3                                       "Compression ratio";

  /******************Instantaneous Fluid Properties***************************/
  SI.Density rho                             "Initial values of density in Kg/m3";
  SI.Temperature T                           "Initial value of temperature in in K";
  SI.Pressure P(displayUnit="kPa")           "Pressure at specified temperature and density conditions";
  SI.SpecificEnthalpy h                      "Enthalpy at specified temperature and density conditions";
  SI.SpecificEnergy u                        "Specific Internal enerygy";
  SI.DynamicViscosity mu=11.668e-6           "Dynamic Viscosity";
  parameter Real Pr=0.8                      "Prandtle number";

  /******************Suction Fluid Properties***************************/

  SI.Pressure P_s                           "Suction Pressure";
  parameter SI.Temperature T_s=297          "Suction Tempersture";
  parameter SI.Density rho_s=9.299           "Suction Density";
  SI.SpecificEnthalpy h_s                   "suction flow specific Enthalpy";

  /******************Discharge Fluid Properties***************************/
  SI.Pressure P_d;

  /******************ThermoPhysical Properties***************************/
  parameter SI.RatioOfSpecificHeatCapacities gamma=1.239       "Specific heat capacity ratio";
  parameter Real R=81.49                                       "Gas Constant";
  Real drhodT=-6.1634e3;
  Real dudT=5.4693e3;

  /****************** Outside Ambient Conditions ***************************/
  parameter SI.DynamicViscosity mu_air=16.82e-6;
  parameter SI.ThermalConductivity k_solid=160   "Shell thermal conductivity";
  parameter SI.Density rho_air=1.2               "Air density";
  parameter SI.Velocity V_air=1                  "Wind velocity";
  parameter Real Pr_air=0.7                      "Outside air Prandtle number";

  /****************** Mass Flow Rates ***************************/
  SI.MassFlowRate mdot_in;
  SI.MassFlowRate mdot_out;
  SI.MassFlowRate mdot;

  /****************** Intermediate Variables ***************************/
  SI.Angle theta_01(displayUnit="Deg")                           "Crank angle in rad";
  SI.Length x_calc                                               "piston position at each time step";
  SI.Length x_2;
  SI.Area A_piston                                               "piston area or cylinder cross sectional area";
  SI.Volume V(displayUnit="cm3")                                 "Volume at each time step inside the cylinder";
  SI.VolumeFlowRate dVdt(displayUnit="cm3/s")                   "Volume derivative w.r.t time";
  SI.HeatFlowRate Qdot;
  SI.Temperature Tamb=295;
  SI.Temperature T_w;
  SI.Velocity  V_s                                              "Sunction velocity";
  SI.Area      A_s;

  /**************************** Valve Inputs **********************************/
  SI.Velocity V_port(start=1)              "Velocity in the port";

  parameter SI.Length d_port=0.0059              "Port diameter";
  parameter SI.Length d_valve=0.007              "Valve diamater in";
  Real C_D=1.17                                  "Drang Coefficient";
  parameter SI.ModulusOfElasticity E=1.93e11     "youngs modulus";
  parameter SI.Length h_valve= 0.0001532         "Valve thickness";
  parameter SI.Length l_valve=0.018              "Total length of the valve";
  parameter SI.Length a_valve=0.0140             "distance from anchor to force";
  parameter SI.Density rho_valve=8000            "Density of the spring steel";
  parameter SI.Length h_disStopper=0.0018        "Discharge stopper height";

  /********************* valve Intermediate variables***************************/

  SI.Length x_tr                                 "transition valve lift";
  SI.Area A_port                                 "Port cross sectional area";
  Real k_valve(unit="N/m")                       "Valve stiffness";
  SI.MomentOfInertia I                           "moment of inertia of the vlalve";
  SI.Mass m_eff                                  "effective mass of the valve";
  SI.Length x_v(max=0.0018)                                    "valve lift";
  SI.Velocity x1_v                                 "First derivative of the valve lift";
  SI.Acceleration x2_v                             "Second derivative of the valve lift";

  /****************** VLE Fluid import ***************************/

TILMedia.VLEFluid_dT vleFluid(
  computeTransportProperties=true,
                              T=T,d=rho,
  final vleFluidType=sim.vleFluidType1)
  annotation (Placement(transformation(extent={{-12,46},{8,66}})));
TILMedia.VLEFluid_dT vleFluid1(
  computeTransportProperties=true,
                               T=T_s,d=rho_s,
  final vleFluidType=sim.vleFluidType1)
  annotation (Placement(transformation(extent={{-14,-2},{6,18}})));

  inner TIL.SystemInformationManager sim(redeclare
      TILMedia.VLEFluidTypes.Refprop_CO2 vleFluidType1)
    annotation (Placement(transformation(extent={{54,66},{74,86}})));

initial equation
  T=T_s;
  rho=rho_s;
  x_v=0;

equation

  // Fluid properties of refrigerant mixture inside the cmpressor
  P=vleFluid.p;
  h=vleFluid.h;

  u=h-P/rho;

  //suction side fluid properties

  P_s=vleFluid1.p;
  h_s=vleFluid1.h;

  // Discharge side Pressure
  P_d=P_s*PR;

  /*************************** Volume Calculation ***************************/
  der(theta_01)= w;
  x_calc = R_01*Modelica.Math.cos(theta_01) + (R_12^2 - R_01^2*
  Modelica.Math.sin(theta_01)^ 2)^0.5;

  x_2 = -x_calc+(R_12 + R_01);
  A_piston = Modelica.Constants.pi*(B/2)^2;
  V = x_2*A_piston+V_dead;
  dVdt=der(V);

  /********************************** Mass Flow *******************************/

  V_s=functions.suction_velocity(P,P_s,T_s,rho_s,gamma,R);
  A_s=Modelica.Constants.pi*d*x_v;
  mdot_in=rho_s*V_s*A_s;
  mdot_out=functions.discharge_valve(P,P_d,T,rho,gamma,R,d);
  mdot=mdot_in-mdot_out;

  /***************************** Heat Transfer ***************************/

  Qdot=functions.Inside_HT(T,rho,T_w,V,dVdt,B,mu,Pr,13.5);
  T_w=functions.outside_HT(Tamb,mu_air,k_solid,rho_air,V_air,D_shell,L_shell,Pr_air);

  /************************* Mass Energy Balance **********************/

  der(rho)=(1/V)*(-rho*dVdt+(mdot_in-mdot_out));
  der(T)=((-rho*h*dVdt)-((u*V+rho*V*drhodT)*der(rho))+(Qdot+mdot_in*h_s-mdot_out*h))/(rho*V*dudT);

  /**********************************  Suction Valve Model ******************/
  x_tr=0.25*(d_port^2/d_valve);
  I=(d_valve*h_valve^3)/12;
  k_valve=(6*E*I)/(a_valve^2*(3*l_valve-a_valve));
  m_eff=(1/3)*rho_valve*l_valve*d_valve*h_valve;
  A_port=3.14*d_port^2/4;
  V_port=mdot_in/(A_port*rho_s);

  //if x_v>h_disStopper then
    //x_v=0.0018;
  //else
    //x_v=x_v;
  //end if;

  if x_v>x_tr or x_v==x_tr then
    x1_v=der(x_v);
    x2_v=der(x1_v);
    m_eff*x2_v+k_valve*x_v=rho_s*(V_port-x1_v)^2*A_port+0.5*C_D*rho_s*(V_port^2-x1_v^2)*A_port;

  else
    x1_v=der(x_v);
    x2_v=der(x1_v);
    m_eff*x2_v+k_valve*x_v=(P_s-P)*A_port+0.5*C_D*rho_s*(V_port^2-x1_v^2)*A_port;
  end if;

annotation (experiment(
      StopTime=0.06,
      __Dymola_NumberOfIntervals=5000,
      __Dymola_Algorithm="Rkfix4"));
end Rec_Compressor_v;
end complete_model;

package blocks
  " Compressor model block. take multiple inputs and outputs"
  model Comp_block
    "This is the coplete block for piston cylinder compressor. It takes three inputs in following order. 
 suction temperature, suction pressure and discharge pressure. output of the block is discharge mass flow rate."
    /*************************Extends********************************/
    extends Modelica.Blocks.Interfaces.MISO(nin=2);
    extends Modelica.Icons.RotationalSensor;
    extends Modelica.Icons.UnderConstruction;


    /****************************Imports*************************************/
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;

    /*********************Compressor Specifications******************************/
    parameter SI.Length R_01=(0.0115)                         "Crank radius in m";
    parameter SI.Length R_12=(0.1055)                         "connecting rod length";
    parameter SI.Length B=(0.02)                              "cylinder bore diamter in m";
    parameter SI.Volume V_dead=(0.00000038)                   "clearance Volume in m3";
    parameter SI.Diameter D_shell=0.05;
    parameter SI.Length L_shell=0.254;
    parameter SI.Diameter d=0.0059                            "Valve diameters";
    parameter SI.AngularVelocity w =   from_rpm(3000);
    parameter Real PR=3                                       "Compression ratio";

    /******************Instantaneous Fluid Properties***************************/
    SI.Density rho                             "Initial values of density in Kg/m3";
    SI.Temperature T                           "Initial value of temperature in in K";
    SI.Pressure P(displayUnit="kPa")           "Pressure at specified temperature and density conditions";
    SI.SpecificEnthalpy h                      "Enthalpy at specified temperature and density conditions";
    SI.SpecificEnergy u1                        "Specific Internal enerygy";
    SI.DynamicViscosity mu=11.668e-6           "Dynamic Viscosity";
    parameter Real Pr=0.8                      "Prandtle number";

    /******************Suction Fluid Properties***************************/

    SI.Pressure P_s=u[2]                           "Suction Pressure";
    SI.Temperature T_s=u[1]          "Suction Tempersture";
    SI.Density rho_s          "Suction Density";
    SI.SpecificEnthalpy h_s                   "suction flow specific Enthalpy";

    /******************Discharge Fluid Properties***************************/
    SI.Pressure P_d;

    /******************ThermoPhysical Properties***************************/
    parameter SI.RatioOfSpecificHeatCapacities gamma=1.239       "Specific heat capacity ratio";
    parameter Real R=81.49                                       "Gas Constant";
    Real drhodT=-6.1634e3                         "Density derivative w.r.t temperature";
    Real dudT=5.4693e3                            "internal energy derivate w.r.t temperature";

    /****************** Outside Ambient Conditions ***************************/
    parameter SI.DynamicViscosity mu_air=16.82e-6;
    parameter SI.ThermalConductivity k_solid=160   "Shell thermal conductivity";
    parameter SI.Density rho_air=1.2               "Air density";
    parameter SI.Velocity V_air=1                  "Wind velocity";
    parameter Real Pr_air=0.7                      "Outside air Prandtle number";

    /****************** Mass Flow Rates ***************************/
    SI.MassFlowRate mdot_in;
    SI.MassFlowRate mdot_out=y;
    SI.MassFlowRate mdot;

    /****************** Intermediate Variables ***************************/
    SI.Angle theta_01(displayUnit="Deg")                           "Crank angle in rad";
    SI.Length x_calc                                               "piston position at each time step";
    SI.Length x_2;
    SI.Area A_piston                                               "piston area or cylinder cross sectional area";
    SI.Volume V(displayUnit="cm3")                                 "Volume at each time step inside the cylinder";
    SI.VolumeFlowRate dVdt(displayUnit="cm3/s")                   "Volume derivative w.r.t time";
    SI.HeatFlowRate Qdot;
    SI.Temperature Tamb=295;
    SI.Temperature T_w;

    /****************** VLE Fluid import ***************************/

  TILMedia.VLEFluid_dT vleFluid(
    computeTransportProperties=true,
                                T=T,d=rho,
    final vleFluidType=sim.vleFluidType1)
    annotation (Placement(transformation(extent={{-12,46},{8,66}})));
  TILMedia.VLEFluid_pT vleFluid1(
    computeTransportProperties=true,
                                 T=T_s,p=P_s,
    final vleFluidType=sim.vleFluidType1)
    annotation (Placement(transformation(extent={{-14,-2},{6,18}})));

    inner TIL.SystemInformationManager sim(redeclare
        TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType1)
      annotation (Placement(transformation(extent={{54,66},{74,86}})));
  initial equation
    T=T_s;
    rho=rho_s;


  equation

    // Fluid properties of refrigerant mixture inside the cmpressor
    P=vleFluid.p;
    h=vleFluid.h;

    u1=h-P/rho;

    //suction side fluid properties

    rho_s=vleFluid1.d;
    h_s=vleFluid1.h;

    // Discharge side Pressure
    P_d=P_s*PR;

    /*************************** Volume Calculation ***************************/
    der(theta_01)= w;
    x_calc = R_01*Modelica.Math.cos(theta_01) + (R_12^2 - R_01^2*
    Modelica.Math.sin(theta_01)^ 2)^0.5;

    x_2 = -x_calc+(R_12 + R_01);
    A_piston = Modelica.Constants.pi*(B/2)^2;
    V = x_2*A_piston+V_dead;
    dVdt=der(V);

    /********************************** Mass Flow *******************************/

    mdot_in=functions.suction_valve(P,P_s,T_s,rho_s,gamma,R,d);
    mdot_out=functions.discharge_valve(P,P_d,T,rho,gamma,R,d);
    mdot=mdot_in-mdot_out;

    /***************************** Heat Transfer ***************************/

    Qdot=functions.Inside_HT(T,rho,T_w,V,dVdt,B,mu,Pr,13.5);
    T_w=functions.outside_HT(Tamb,mu_air,k_solid,rho_air,V_air,D_shell,L_shell,Pr_air);

    /************************* Mass Energy Balance **********************/

    der(rho)=(1/V)*(-rho*dVdt+(mdot_in-mdot_out));
    der(T)=((-rho*h*dVdt)-((u1*V+rho*V*drhodT)*der(rho))+(Qdot+mdot_in*h_s-mdot_out*h))/(rho*V*dudT);


  annotation (experiment(
        StopTime=0.06,
        __Dymola_NumberOfIntervals=1000,
        __Dymola_Algorithm="Cerk45"));
  end Comp_block;

package models "Valve and volume model"
  model valve_model "This model is based on 1-D mass spring vibration model. It 
  simulates time dependent valve opening in response to the pressure difference. 
  Also, it divides valve opening into two regions i.e. pressure dominant and 
  mass flux dominant region."

    /**************************** Imports **********************************/
    import SI = Modelica.SIunits;

    /**************************** Inputs **********************************/
    parameter SI.Pressure P_high=400            "High pressure side pressure";
    parameter SI.Pressure P_low(start=300)             "Low pressure side pressure";
    parameter SI.Density rho_port(start=10)            "Desnity of the fluid at port conditions";
    parameter SI.Velocity V_port(start=1)              "Velocity in the port";

    parameter SI.Length d_port=0.0059              "Port diameter";
    parameter SI.Length d_valve=0.007              "Valve diamater in";
    Real C_D=1.17                                  "Drang Coefficient";
    parameter SI.ModulusOfElasticity E=1.93e11     "youngs modulus";
    parameter SI.Length h_valve= 0.0001532         "Valve thickness";
    parameter SI.Length l_valve=0.018              "Total length of the valve";
    parameter SI.Length a_valve=0.0140             "distance from anchor to force";
    parameter SI.Density rho_valve=8000            "Density of the spring steel";
    parameter SI.Length h_disStopper=0.0018        "Discharge stopper height";

    /********************* Intermediate variables***************************/

    SI.Length x_tr                                 "transition valve lift";
    SI.Area A_port                                 "Port cross sectional area";
    Real k_valve(unit="N/m")                       "Valve stiffness";
    SI.MomentOfInertia I                           "moment of inertia of the vlalve";
    SI.Mass m_eff                                  "effective mass of the valve";
    SI.Length x_v                                    "valve lift";
    SI.Velocity x1_v                                 "First derivative of the valve lift";
    SI.Acceleration x2_v                             "Second derivative of the valve lift";

  initial equation
    x_v=0;

  equation
    x_tr=0.25*(d_port^2/d_valve);
    I=(d_valve*h_valve^3)/12;
    k_valve=(6*E*I)/(a_valve^2*(3*l_valve-a_valve));
    m_eff=(1/3)*rho_valve*l_valve*d_valve*h_valve;
    A_port=3.14*d_port^2/4;
    if x_v>x_tr or x_v==x_tr then
      x1_v=der(x_v);
      x2_v=der(x1_v);
      m_eff*x2_v+k_valve*x_v=rho_port*(V_port-x1_v)^2*A_port+0.5*C_D*rho_port*V_port^2*A_port;
    else
      x1_v=der(x_v);
      x2_v=der(x1_v);
      m_eff*x2_v+k_valve*x_v=(P_high-P_low)*A_port+0.5*C_D*rho_port*V_port^2*A_port;
    end if;

    annotation (experiment(StopTime=0.3, __Dymola_Algorithm="Esdirk45a"));
  end valve_model;

  model Volume "This model calculates the volume and volume derivative inside the 
  cylinder at each time instance"

    /************************** Imports ****************************/
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;

    /************************** Inputs ****************************/

    input SI.Volume V_disp;
    input SI.AngularVelocity w     "Angular Velocity";
    input SI.Volume V_dead         "clearance Volume in m3";

    /****************** Intermediate variables *******************/
    SI.Angle theta_01              " Crank angle";

    /************************** outputs ****************************/
    output SI.Volume V             "Instantaneous volume inside cylinder";

    output SI.VolumeFlowRate dVdt(displayUnit="cm3/s")    "Volume derivative w.r.t time";

  equation
    der(theta_01)= w;

  V = V_dead+V_disp/2*(1-cos(theta_01));
  dVdt=der(V);
  end Volume;
end models;
end blocks;

package testers
  model Comp_block_tester

      Compressor_model.blocks.Comp_block comp_block
        annotation (Placement(transformation(extent={{4,0},{24,20}})));
    Modelica.Blocks.Sources.Constant const1(k=508)
      annotation (Placement(transformation(extent={{-72,-28},{-52,-8}})));
      inner TIL.SystemInformationManager sim(redeclare
          TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType1)
        annotation (Placement(transformation(extent={{44,60},{64,80}})));
    Modelica.Blocks.Sources.Constant const2(k=297)
      annotation (Placement(transformation(extent={{-76,32},{-56,52}})));
  equation
      connect(const2.y, comp_block.u[1]) annotation (Line(points={{-55,42},{-26,
              42},{-26,10},{2,10}}, color={0,0,127}));
      connect(const1.y, comp_block.u[2]) annotation (Line(points={{-51,-18},{
              -24,-18},{-24,10},{2,10}}, color={0,0,127}));
      annotation ();
  end Comp_block_tester;
end testers;

package in_progress
  " these models,blocks or testers do not work properly and are under development phase"
  model CV1_block_1
      "This model calculates the Volume of the Control volume and change in control volume w.r.t time. Volume model is changed in this code. instead of entering the dimensions of the ocmpressor, it takes only maximum displaced volume as input"

    extends Modelica.Icons.UnderConstruction;
  import SI = Modelica.SIunits;
  import Modelica.SIunits.Conversions.*;
  parameter SI.Volume V_disp=7.6057e-6  "totoal displaced volume of the compressor";
  // Compressor dimensions
     parameter SI.Length B=(0.02)                              "cylinder bore diamter in m";
     parameter SI.Volume V_dead=(0.00000038)                   "clearance Volume in m3";
     parameter SI.Length d=0.0059                              "Valve diameters";

  // Compressor RPM

  parameter SI.AngularVelocity w =   from_rpm(3000);
  parameter Real PR=3                              " Pressure ratio";

  //**********************   Connectors ******************************

    public
           TIL.Connectors.VLEFluidPort portA(final vleFluidType(
        fixedMixingRatio=false,
        nc_propertyCalculation=2,
        vleFluidNames={"TILMedia.R134A"},
        mixingRatio_propertyCalculation={1}))
      "Suction port"
      annotation (Placement(transformation(extent={{-10,-90},{10,-70}}, rotation=
              0)));
  TIL.Connectors.VLEFluidPort portB(final vleFluidType(
        fixedMixingRatio=false,
        nc_propertyCalculation=2,
        vleFluidNames={"TILMedia.R134A"},
        mixingRatio_propertyCalculation={1}))
      "Discharge port"
      annotation (Placement(transformation(extent={{-10,70},{10,90}}, rotation=0)));
  // Fluid Properties
  SI.Density rho                              "Initial values of density in Kg/m3";
  SI.Temperature T                            "Initial value of temperature in in K";
  SI.Pressure P(displayUnit="kPa")           "Pressure at specified temperature and density conditions";
  SI.SpecificEnthalpy h                      "Enthalpy at specified temperature and density conditions";
  SI.SpecificEnergy u2                        "Specific Internal enerygy";
  SI.DynamicViscosity mu=11.668e-6;
  Real Pr=0.8                    "Prandtle number";

  //suction side fluid properties

  SI.Pressure P_s=portA.p;
  SI.Temperature T_s;
  SI.Density rho_s;
  SI.SpecificEnthalpy h_s=portA.h_outflow                      "suction flow specific Enthalpy";

  //discharge side fluid properties

  SI.Pressure P_d;

  Real gamma=1.239;
  Real R=81.49;

  // outside air properties
  SI.DynamicViscosity mu_air=16.82e-6;
  SI.ThermalConductivity k_solid=160;
  SI.Density rho_air=1.2;
  SI.Velocity V_air=1;
  Real Pr_air=0.7;

  // Mass flow rates through valves
  SI.MassFlowRate mdot_in;
  SI.MassFlowRate mdot_out;
  SI.MassFlowRate mdot;
  SI.Mass mdot_tot_in;
  SI.Mass mdot_tot_out;

    // Volume calculation variables
  SI.Angle theta_01(displayUnit="Deg")                           "Crank angle in rad";
  // SI.Length x_calc                                               "piston position at each time step";
  // SI.Length x_2;
  // SI.Area A_piston                                               "piston area or cylinder cross sectional area";
  SI.Volume V(displayUnit="cm3")                                 "Volume at each time step inside the cylinder";
  SI.VolumeFlowRate dVdt(displayUnit="cm3/s")                   "Volume derivative w.r.t time";
  SI.HeatFlowRate Qdot;
  //SI.HeatFlowRate Qdot_out;
  SI.Temperature Tamb=295;
  SI.Temperature T_w;

  SI.Power Wdot;
  SI.Diameter D_shell=0.05;
  SI.Length L_shell=0.254;

  TILMedia.VLEFluid_dT vleFluid(
    computeTransportProperties=true,
                                T=T,d=rho,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType)
    annotation (Placement(transformation(extent={{-12,46},{8,66}})));
  TILMedia.VLEFluid_dT vleFluid1(
    computeTransportProperties=true,
                                 T=T_s,d=rho_s,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType)
    annotation (Placement(transformation(extent={{-14,-2},{6,18}})));

    TILMedia.VLEFluid_ph vleFluid2(
    computeTransportProperties=true,
                                 p=P_s,h=h_s,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType)
      annotation (Placement(transformation(extent={{-14,-42},{6,-22}})));
  initial equation
  T=297;
  rho=23.7;

  equation

  // Fluid properties of refrigerant mixture inside the cmpressor
  P=vleFluid.p;
  h=vleFluid.h;

  u2=h-P/rho;

  //suction side fluid properties

  T_s=vleFluid2.T;
  rho_s=vleFluid2.d;

  // Discharge side Pressure
  P_d=P_s*PR;

  //This section of the code is for calculation of Volume and derivative of Volume w.r.t
  der(theta_01)= w;

  V = V_dead+V_disp/2*(1-cos(theta_01));
  dVdt=der(V);

  // Mass flow calculations

  mdot_in=suction_valve(P,P_s,T_s,rho_s,gamma,R,d);
  mdot_out=discharge_valve(P,P_d,T,rho,gamma,R,d);

  mdot=mdot_in-mdot_out;
  Qdot=Inside_HT(T,rho,T_w,V,dVdt,B,mu,Pr,13.5);
  T_w=outside_HT(Tamb,mu_air,k_solid,rho_air,V_air,D_shell,L_shell,Pr_air);

  der(rho)=(1/V)*(-rho*dVdt+(mdot_in-mdot_out));
  der(T)=((-rho*h*dVdt)-((u2*V+rho*V*(-6.1634e3))*der(rho))+(Qdot+mdot_in*h_s-mdot_out*h))/(rho*V*5.4693e3);

  der(mdot_tot_in)=mdot_in;
  der(mdot_tot_out)=mdot_out;

  Wdot=P*der(V);
  portA.m_flow=mdot_in;
  portB.m_flow=mdot_out;
  portB.h_outflow=h;
  portB.p=P;
  annotation (experiment(
        StopTime=0.02,
        __Dymola_NumberOfIntervals=2000,
        __Dymola_Algorithm="Cerk45"));
  end CV1_block_1;

  model CV1_block_1_2 "This model calculates the Volume of the Control volume and change in control volume 
  w.r.t time. Volume model is changed in this code. 

  instead of entering the dimensions of the ocmpressor, it takes only maximum displaced volume as input"

    extends Modelica.Icons.UnderConstruction;
  import SI = Modelica.SIunits;
  import Modelica.SIunits.Conversions.*;
  parameter SI.Volume V_disp=7.6057e-6  "totoal displaced volume of the compressor";
  // Compressor dimensions
     parameter SI.Length B=(0.02)                              "cylinder bore diamter in m";
     parameter SI.Volume V_dead=(0.00000038)                   "clearance Volume in m3";
     parameter SI.Length d=0.0059                              "Valve diameters";

  // Compressor RPM

  parameter SI.AngularVelocity w =   from_rpm(3000);
  parameter Real PR=3                              " Pressure ratio";

  //**********************   Connectors ******************************

    public
           TIL.Connectors.VLEFluidPort portA(final vleFluidType = sim.vleFluidType1)
      "Suction port"
      annotation (Placement(transformation(extent={{-10,-90},{10,-70}}, rotation=
              0)));
  TIL.Connectors.VLEFluidPort portB(final vleFluidType = sim.vleFluidType1)
      "Discharge port"
      annotation (Placement(transformation(extent={{-10,70},{10,90}}, rotation=0)));
  // Fluid Properties
  SI.Density rho(start=23.7)                              "Initial values of density in Kg/m3";
  SI.Temperature T(start=297)                            "Initial value of temperature in in K";
  SI.Pressure P(displayUnit="kPa",start=500000)           "Pressure at specified temperature and density conditions";
  SI.SpecificEnthalpy h                            "Enthalpy at specified temperature and density conditions";
  SI.SpecificEnergy u2                        "Specific Internal enerygy";
  SI.DynamicViscosity mu=11.668e-6;
  Real Pr=0.8                    "Prandtle number";

  //suction side fluid properties

  SI.Pressure P_s;
  SI.Temperature T_s;
  SI.Density rho_s;
  SI.SpecificEnthalpy h_s                      "suction flow specific Enthalpy";

  //discharge side fluid properties

  SI.Pressure P_d(min=1);

  Real gamma=1.239;
  Real R=81.49;

  // outside air properties
  SI.DynamicViscosity mu_air=16.82e-6;
  SI.ThermalConductivity k_solid=160;
  SI.Density rho_air=1.2;
  SI.Velocity V_air=1;
  Real Pr_air=0.7;

  // Mass flow rates through valves
  SI.MassFlowRate mdot_in;
  SI.MassFlowRate mdot_out;
  SI.MassFlowRate mdot;

    // Volume calculation variables
  SI.Angle theta_01(displayUnit="Deg")                           "Crank angle in rad";
  // SI.Length x_calc                                               "piston position at each time step";
  // SI.Length x_2;
  // SI.Area A_piston                                               "piston area or cylinder cross sectional area";
  SI.Volume V(displayUnit="cm3")                                 "Volume at each time step inside the cylinder";
  SI.VolumeFlowRate dVdt(displayUnit="cm3/s")                   "Volume derivative w.r.t time";
  SI.HeatFlowRate Qdot;
  //SI.HeatFlowRate Qdot_out;
  SI.Temperature Tamb=295;
  SI.Temperature T_w;

  SI.Diameter D_shell=0.05;
  SI.Length L_shell=0.254;

  // ******************************** Fluid Properties *********************
  TILMedia.VLEFluid_dT vleFluid(
    computeTransportProperties=true,
                                T=T,d=rho,
      final vleFluidType=sim.vleFluidType1)
    annotation (Placement(transformation(extent={{-12,46},{8,66}})));

    TILMedia.VLEFluid_ph vleFluid2(
    computeTransportProperties=true,
                                 p=portA.p,h=inStream(portA.h_outflow),
      final vleFluidType= sim.vleFluidType1)
      annotation (Placement(transformation(extent={{-14,-42},{6,-22}})));
    inner TIL.SystemInformationManager sim(redeclare
        TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType1)
      annotation (Placement(transformation(extent={{50,64},{70,84}})));
  equation

  // Fluid properties of refrigerant mixture inside the cmpressor
  P=vleFluid.p;
  h=vleFluid.h;

  u2=h-P/rho;

  //suction side fluid properties
  P_s=-portA.p;
  h_s=portA.h_outflow;
  T_s=vleFluid2.T;
  rho_s=vleFluid2.d;

  // Discharge side Pressure
  P_d=P_s*PR;
  //P_d=portB.p;

  //This section of the code is for calculation of Volume and derivative of Volume w.r.t
  der(theta_01)= w;

  V = V_dead+V_disp/2*(1-cos(theta_01));
  dVdt=der(V);

  // Mass flow calculations

  mdot_in=suction_valve(P,P_s,T_s,rho_s,gamma,R,d);
  mdot_out=discharge_valve(P,P_d,T,rho,gamma,R,d);

  mdot=mdot_in-mdot_out;
  Qdot=Inside_HT(T,rho,T_w,V,dVdt,B,mu,Pr,13.5);
  T_w=outside_HT(Tamb,mu_air,k_solid,rho_air,V_air,D_shell,L_shell,Pr_air);

  der(rho)=(1/V)*(-rho*dVdt+(mdot_in-mdot_out));
  der(T)=((-rho*h*dVdt)-((u2*V+rho*V*(-6.1634e3))*der(rho))+(Qdot+mdot_in*h_s-mdot_out*h))/(rho*V*5.4693e3);

  portA.m_flow=mdot_in;
  portB.m_flow=mdot_out;
  portB.h_outflow=h_s;
  portA.h_limit = -1e6; //no cell volume!
  portB.p=P_d;
    portB.h_limit = -1e6;
    portA.xi_outflow = inStream(portB.xi_outflow);
    portB.xi_outflow = inStream(portA.xi_outflow);
  annotation (experiment(
        StopTime=0.02,
        __Dymola_NumberOfIntervals=2000,
        __Dymola_Algorithm="Dassl"));
  end CV1_block_1_2;

  model test1

    inner TIL.SystemInformationManager sim(
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType1,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType2,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType3,
      redeclare TILMedia.GasTypes.VDI4670_MoistAir gasType1)
                                                     annotation (Placement(
          transformation(extent={{70,70},{90,90}},   rotation=0)));
    TIL.VLEFluidComponents.PressureStateElements.PressureState pressureState_hp(
        pressureStateID=1, pInitial=11e5)
      annotation (Placement(transformation(extent={{-36,54},{-24,66}},
                                                                    rotation=0)));
    TIL.VLEFluidComponents.Separators.Separator separator(pressureStateID=1, V=
          0.5e-3) annotation (Placement(transformation(extent={{-64,46},{-76,66}},
            rotation=0)));
    TIL.HeatExchangers.MPET.MoistAirVLEFluid.CrossFlowHX condenser(
      pressureStateID=1,
      pVLEFluidStart=pressureState_hp.pInitial,
      redeclare model WallHeatConductionModel =
          TIL.HeatExchangers.MPET.TransportPhenomena.WallHeatTransfer.GeometryBasedConduction,
      initVLEFluid="linearEnthalpyDistribution",
      m_flowVLEFluidStart=0.02,
      nCellsPerPass=3,
      redeclare TIL.HeatExchangers.MPET.Geometry.Example hxGeometry(
        hxHeight=0.34,
        nPasses=3,
        nTubesPerPass={16,11,6}),
      hInitialVLEFluid_Cell1=300e3,
      hInitialVLEFluid_CellN=250e3,
      TInitialWall=273.15 + 45,
      redeclare model FinSideHeatTransferModel =
          TIL.HeatExchangers.MPET.TransportPhenomena.FinSideHeatTransfer.Chang,
      redeclare model TubeSideHeatTransferModel =
          TIL.HeatExchangers.MPET.TransportPhenomena.TubeSideHeatTransfer.ConstantAlpha
            (
           constantAlpha=2000),
      redeclare model FinSidePressureDropModel =
          TIL.HeatExchangers.MPET.TransportPhenomena.FinSidePressureDrop.KimBullard,
      wallCellStateType="state south",
      fixedPressureDropInitialMoistAir=true) annotation (Placement(transformation(
            extent={{14,46},{-14,74}},rotation=0)));

    TIL.GasComponents.Boundaries.BoundaryUnderdetermined moistAirSink_condenser
      annotation (Placement(transformation(
          origin={0,36},
          extent={{-4,-10},{4,10}},
          rotation=90)));
    TIL.GasComponents.Boundaries.BoundaryOverdetermined moistAirSource_condenser(
      boundaryType="p, m_flow",
      TFixed=273.15 + 30,
      m_flowFixed=-0.2/38*33,
      streamVariablesInputTypeConcentration="phi",
      phiFixed=40) annotation (Placement(transformation(
          origin={0,84},
          extent={{-4,-10},{4,10}},
          rotation=90)));
    TIL.VLEFluidComponents.Valves.OrificeValve valve(use_effectiveFlowAreaInput=
         false, effectiveFlowAreaFixed=6.82320e-7)  annotation (Placement(
          transformation(
          origin={-70,0},
          extent={{-8,-4},{8,4}},
          rotation=270)));
    TIL.VLEFluidComponents.PressureStateElements.PressureState pressureState_lp(
        pressureStateID=2, pInitial=250000) annotation (Placement(transformation(
            extent={{-36,-36},{-24,-24}},
                                      rotation=0)));
    TIL.VLEFluidComponents.Tubes.Tube tube(
      pressureStateID=2,
      redeclare TIL.VLEFluidComponents.Tubes.Geometry.TubeGeometry tubeGeometry(
          crossSectionType=TIL.Internals.CrossSectionType.Circular),
      redeclare model TubeSideHeatTransferModel =
          TIL.VLEFluidComponents.Tubes.TransportPhenomena.HeatTransfer.ConstantAlpha
            (
           constantAlpha=0),
      redeclare model WallHeatConductionModel =
          TIL.VLEFluidComponents.Tubes.TransportPhenomena.WallHeatTransfer.GeometryBasedConduction,
      redeclare model WallMaterial = TILMedia.SolidTypes.TILMedia_Steel,
      nCells=3,
      redeclare model PressureDropModel =
          TIL.VLEFluidComponents.Tubes.TransportPhenomena.PressureDrop.ConstantZeta
            (
           zeta=5, massFlowLimit=1e-12))
      annotation (Placement(transformation(
          extent={{-8,-2},{8,2}},
          rotation=0,
          origin={-50,-30})));

    TIL.HeatExchangers.Plate.VLEFluidLiquid.ParallelFlowHX parallelFlowHX(
      pressureStateID=2,
      nCells=5,
      redeclare model HeatTransferModel_a =
          TIL.HeatExchangers.Plate.TransportPhenomena.HeatTransfer.ConstantAlpha
            (
           constantAlpha=2000),
      redeclare model PressureDropModel_a =
          TIL.HeatExchangers.Plate.TransportPhenomena.PressureDrop.ZeroPressureDrop,
      cellOrientation_a="B",
      redeclare model HeatTransferModel_b =
          TIL.HeatExchangers.Plate.TransportPhenomena.HeatTransfer.ConstantAlpha
            (
           constantAlpha=300),
      redeclare model PressureDropModel_b =
          TIL.HeatExchangers.Plate.TransportPhenomena.PressureDrop.ZeroPressureDrop,
      redeclare model WallMaterial = TILMedia.SolidTypes.TILMedia_Steel,
      redeclare model WallHeatConductionModel =
          TIL.HeatExchangers.Plate.TransportPhenomena.WallHeatTransfer.GeometryBasedConduction,
      pressureDropInitial_a=0,
      redeclare TIL.HeatExchangers.Plate.Geometry.Example hxGeometry,
      fixedPressureDropInitial_a=true,
      initVLEFluid_a="linearEnthalpyDistribution",
      initLiquid_b="linearTemperatureDistribution",
      hInitialVLEFluid_a_Cell1=410e3,
      hInitialVLEFluid_a_CellN=380e3,
      TInitialLiquid_b_Cell1=278.15,
      TInitialLiquid_b_CellN=283.15) annotation (Placement(transformation(
          extent={{-14,-14},{14,14}},
          rotation=180,
          origin={0,-38})));

    TIL.LiquidComponents.Boundaries.Boundary boundary(boundaryType="p")
      annotation (Placement(transformation(extent={{-32,-60},{-24,-40}})));
    TIL.LiquidComponents.Boundaries.Boundary boundary1(boundaryType="p")
      annotation (Placement(transformation(extent={{52,-60},{60,-40}})));
    TIL.LiquidComponents.Pumps.SimplePump simplePump(V_flowFixed=1/100/60)
      annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={34,-50})));
    TIL.VLEFluidComponents.Sensors.StatePoint statePoint
      annotation (Placement(transformation(extent={{46,64},{54,72}})));
    TIL.VLEFluidComponents.Sensors.StatePoint statePoint1(stateViewerIndex=1)
      annotation (Placement(transformation(extent={{-84,24},{-76,32}})));
    TIL.VLEFluidComponents.Sensors.StatePoint statePoint2(stateViewerIndex=2)
      annotation (Placement(transformation(extent={{-84,-26},{-76,-18}})));
    TIL.VLEFluidComponents.Sensors.StatePoint statePoint3(stateViewerIndex=3)
      annotation (Placement(transformation(extent={{46,-26},{54,-18}})));
      Compressor_model.in_progress.CV1_block_1_2 cV1_block_1_2_1(
        V_disp=150e-6,
        d=0.59,
        w=158) annotation (Placement(transformation(extent={{30,10},{50,30}})));
  equation
    connect(condenser.portB_gas,moistAirSink_condenser. port) annotation (Line(
        points={{0,46},{0,36}},
        color={255,153,0},
        thickness=0.5));
    connect(separator.portInlet,pressureState_hp. portB) annotation (Line(
        points={{-65,60},{-36,60}},
        color={153,204,0},
        thickness=0.5));
    connect(pressureState_hp.portA,condenser. portB_vle) annotation (Line(
        points={{-24,60},{-14,60}},
        color={153,204,0},
        thickness=0.5));
    connect(separator.portLiquid,valve. portA) annotation (Line(
        points={{-70,46},{-70,8}},
        color={153,204,0},
        thickness=0.5));
    connect(valve.portB,tube. portA) annotation (Line(
        points={{-70,-8},{-70,-30},{-58,-30}},
        color={153,204,0},
        thickness=0.5));
    connect(pressureState_lp.portB,tube. portB) annotation (Line(
        points={{-36,-30},{-42,-30}},
        color={153,204,0},
        thickness=0.5));
    connect(pressureState_lp.portA,parallelFlowHX. portB_a) annotation (Line(
        points={{-24,-30},{-14,-30}},
        color={153,204,0},
        thickness=0.5));
    connect(parallelFlowHX.portA_b,simplePump. portA) annotation (Line(
        points={{14,-46},{20,-46},{20,-50},{26,-50}},
        color={0,170,238},
        thickness=0.5));
    connect(simplePump.portB,boundary1. port) annotation (Line(
        points={{42,-50},{56,-50}},
        color={0,170,238},
        thickness=0.5));
    connect(parallelFlowHX.portB_b,boundary. port) annotation (Line(
        points={{-14,-46},{-20,-46},{-20,-50},{-28,-50}},
        color={0,170,238},
        thickness=0.5));
    connect(moistAirSource_condenser.port,condenser. portA_gas) annotation (Line(
        points={{0,84},{0,74}},
        color={255,153,0},
        thickness=0.5));
    connect(valve.portA,statePoint1. sensorPort) annotation (Line(
        points={{-70,8},{-70,20},{-80,20},{-80,24}},
        color={153,204,0},
        thickness=0.5));
    connect(tube.portA,statePoint2. sensorPort) annotation (Line(
        points={{-58,-30},{-80,-30},{-80,-26}},
        color={153,204,0},
        thickness=0.5));
    connect(parallelFlowHX.portA_a, cV1_block_1_2_1.portA) annotation (Line(
        points={{14,-30},{40,-30},{40,12}},
        color={153,204,0},
        thickness=0.5));
    connect(statePoint3.sensorPort, cV1_block_1_2_1.portA) annotation (Line(
        points={{50,-26},{40,-26},{40,12}},
        color={153,204,0},
        thickness=0.5));
    connect(condenser.portA_vle, cV1_block_1_2_1.portB) annotation (Line(
        points={{14,60},{40,60},{40,28}},
        color={153,204,0},
        thickness=0.5));
    connect(statePoint.sensorPort, cV1_block_1_2_1.portB) annotation (Line(
        points={{50,64},{46,64},{46,56},{40,56},{40,28}},
        color={153,204,0},
        thickness=0.5));
    annotation ();
  end test1;

  model simple_test

    inner TIL.SystemInformationManager
                                   sim(
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType1,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType2,
      redeclare TILMedia.VLEFluidTypes.TILMedia_R134a vleFluidType3)
                             annotation (Placement(transformation(extent={{60,60},
              {80,80}}, rotation=0)));
    TIL.VLEFluidComponents.Boundaries.BoundaryOverdetermined
      boundaryOverdetermined(
      streamVariablesInputType="h",
      hFixed=415000,
      boundaryType="p, m_flow",
      pFixed(displayUnit="bar") = 500000,
      m_flowFixed=0.2)
      annotation (Placement(transformation(extent={{-60,10},{-52,30}})));
    TIL.VLEFluidComponents.Boundaries.BoundaryUnderdetermined
      boundaryUnderdetermined
      annotation (Placement(transformation(extent={{44,14},{52,34}})));
      in_progress.CV1_block_1_2 cV1_block_1_2_1 annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-8,26})));
  equation
    connect(boundaryOverdetermined.port, cV1_block_1_2_1.portA) annotation (
        Line(
        points={{-56,20},{-38,20},{-38,26},{-16,26}},
        color={153,204,0},
        thickness=0.5));
    connect(cV1_block_1_2_1.portB, boundaryUnderdetermined.port) annotation (
        Line(
        points={{0,26},{24,26},{24,24},{48,24}},
        color={153,204,0},
        thickness=0.5));
    annotation (experiment(__Dymola_Algorithm="Esdirk45a"));
  end simple_test;
end in_progress;
  annotation (uses(TIL(version="3.7.0"), TILMedia(version="3.7.0"),
      Modelica(version="3.2.2")),                                    experiment(
      StopTime=0.03,
      __Dymola_NumberOfIntervals=2000,
      __Dymola_Algorithm="Cerk45"));
end Compressor_model;
