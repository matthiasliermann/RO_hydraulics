within ;
package RO_hydraulics "package for modeling reverse osmosis membrane process"
  model RO_system "Example system for RO process"

  //extends OpenHydraulics.Interfaces.PartialFluidComponent;

    extends OpenHydraulics.Interfaces.PartialFluidCircuit(redeclare
       OpenHydraulics.Fluids.GenericOilSimple oil);
    RO rO1 annotation(Placement(visible = true, transformation(origin = {-22, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenHydraulics.Basic.VarPressureSource source annotation(Placement(visible = true, transformation(origin = {-58, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Cosine cosine1(amplitude = 35e5, freqHz = 2, offset = 36e5)  annotation(Placement(visible = true, transformation(origin = {-118, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenHydraulics.Basic.VarVolumeSource source2 annotation(Placement(visible = true, transformation(origin = {-22, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Cosine cosine2(amplitude = 0e-5, freqHz = 2, offset = -20e-5) annotation(Placement(visible = true, transformation(origin = {-68, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(cosine2.y, source2.control) annotation(Line(points={{-57,-28},{-32,
            -28},{-32,-28},{-32,-28}},                                                                             color = {0, 0, 127}));
    connect(rO1.port1, source2.port) annotation(Line(points={{-22,8.2},{-22,8.2},
            {-22,-18},{-22,-18}},                                                                           color = {255, 0, 0}));
    connect(cosine1.y, source.control) annotation(Line(points={{-107,8},{-68,8},
            {-68,8},{-68,8}},                                                                              color = {0, 0, 127}));
    connect(source.port, rO1.port) annotation(Line(points = {{-58, 18}, {-32, 18}}, color = {255, 0, 0}));
    annotation ();
  end RO_system;

  model RO "Series of segments of one RO module modeled as one element"

  import OpenHydraulics.Interfaces;
   extends OpenHydraulics.Interfaces.PartialFluidComponent;
   Real SpecificEnergy(unit="kWh/m^3");
   RO_first rO_first1(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle rO_middle1(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin = {-30, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle rO_middle2(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin = {10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle rO_middle3(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin = {50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_last rO_last1 annotation(Placement(visible = true, transformation(origin = {90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenHydraulics.Interfaces.FluidPort port annotation(Placement(visible = true, transformation(origin = {-84, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenHydraulics.Interfaces.FluidPort port1 annotation(Placement(visible = true, transformation(origin = {90, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(rO_first1.pin_n4, rO_middle1.pin_n) annotation(Line(points = {{-70, 0}, {-70, -10}, {-48, -10}, {-48, 10}, {-40, 10}}, color = {0, 0, 255}));
    connect(rO_first1.pin_n5, rO_middle1.pin_n1) annotation(Line(points = {{-60, 10}, {-54, 10}, {-54, 28}, {-30, 28}, {-30, 20}}, color = {0, 0, 255}));
    connect(rO_middle3.pin_n4, rO_last1.pin_n) annotation(Line(points = {{50, 0}, {50, 0}, {50, -8}, {74, -8}, {74, 10}, {80, 10}, {80, 10}}, color = {0, 0, 255}));
    connect(rO_middle3.pin_n5, rO_last1.pin_n1) annotation(Line(points={{59.8,10},
            {68,10},{68,28},{90,28},{90,20},{90,20},{90,20}},                                                                                    color = {0, 0, 255}));
    connect(rO_middle2.pin_n5, rO_middle3.pin_n1) annotation(Line(points = {{19.8, 10}, {28, 10}, {28, 28}, {40, 28}, {50, 28}, {50, 20}}, color = {0, 0, 255}));
    connect(rO_middle2.pin_n4, rO_middle3.pin_n) annotation(Line(points = {{10, 0}, {10, 0}, {10, -10}, {32, -10}, {32, 10}, {40, 10}}, color = {0, 0, 255}));
    connect(rO_middle1.pin_n5, rO_middle2.pin_n1) annotation(Line(points = {{-20.2, 10}, {-10, 10}, {-10, 28}, {10, 28}, {10, 20}}, color = {0, 0, 255}));
    connect(rO_middle1.pin_n4, rO_middle2.pin_n) annotation(Line(points = {{-30, 0}, {-30, -10}, {-6, -10}, {-6, 10}, {0, 10}}, color = {0, 0, 255}));
    port.m_flow / oil.density(port.p) = rO_first1.Q_fin;
    port1.p = rO_last1.p_fout;
  //rO_last1.Q_fout=50e-5;
    port1.m_flow / oil.density(port1.p) = -rO_last1.Q_fout;
  //rO_first1.p_fin=25;
    port.p = rO_first1.p_fin;
    SpecificEnergy=rO_first1.Q_fin*rO_first1.p_fout/(rO_last1.Q_pout)/(3.6e6);

    annotation (                                                                Icon(graphics={  Polygon(lineColor = {0, 0, 255}, fillColor = {7, 0, 255},
              fillPattern =                                                                                                                                              FillPattern.Solid, points = {{-100, 100}, {100, 100}, {100, -100}, {100, -100}, {-100, 100}}), Line(points = {{-100, -100}, {-100, 100}, {100, 100}, {100, -100}, {-100, -100}})}), Diagram);
  end RO;

  model RO_middle "middle element of RO membrane"
  protected
    parameter Modelica.SIunits.MolarMass M_Na = 22.9898e-3 "Molar mass of sodium [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_Cl = 35.4527e-3 "Molar mass of chlorine [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_s = M_Na + M_Cl "Molar mass of sodium chloride [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_w = 18.0153e-3 "Molar mass of water [kg/mol]";
    parameter Modelica.SIunits.DiffusionCoefficient D = 1.6e-9 "Salt diffusion coefficient [m2/s]";
    parameter Modelica.SIunits.Temperature T = 273 + 25;
    parameter Modelica.SIunits.MolarHeatCapacity gamma = 8.3086 "Molar gas constant R [J/mol.K]";
    // RO module specific parameters
  public
    parameter Real k_f = 9.6 "friction factor, find experimentally";
    parameter Modelica.SIunits.Height h_b = 0.7e-3 "height of...";
    parameter Modelica.SIunits.Length dx = 0.88 / 5 "Length of element";
    parameter Modelica.SIunits.Length W = 1.43 "width of partial element";
    parameter Integer n_leaves = 3 "number of leaves";
    parameter Modelica.SIunits.Length L_mix = 0.006 "characteristic length of mixing net (spacer)";
    parameter Real K = 0.5 "mixing coefficient of net";
    parameter Real A_w(unit="m/(Pa.s)") = 2.08e-12 "RO membrange solvent permeability coefficient";
    Real B_s(unit="m/s") = 1.11e-7 "permeability coefficient of solute across membrane";

    Modelica.SIunits.MolarDensity c_fin "feed inlet molar salt concentration";
    Modelica.SIunits.Pressure p_fin "feed inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_fin "feed inlet volume flowrate";
    Modelica.SIunits.Density rho_fin "feed inlet density";
    Modelica.SIunits.DynamicViscosity mu_fin "feed inlet dynamic viscosity";

    Modelica.SIunits.Pressure p_fout "feed outlet pressure";
    Modelica.SIunits.MolarDensity c_fout "feed outlet molar salt concentration";
    Modelica.SIunits.VolumeFlowRate Q_fout "feed outlet volume flowrate";
    Modelica.SIunits.Density rho_fout "feed outlet density";

    Modelica.SIunits.MolarDensity c_pin "permeate inlet molar salt concentration";
    Modelica.SIunits.Pressure p_pin "permeate inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pin "permeate inlet volume flowrate";
    Modelica.SIunits.Density rho_pin "permeate inlet density";

    Modelica.SIunits.MolarDensity c_pout "permeate outlet molar salt concentration";
    Modelica.SIunits.Pressure p_pout "permeate outlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pout "permeate outlet volume flow rate";
    Modelica.SIunits.Density rho_pout "permeate outlet density";

    Modelica.SIunits.MolarDensity c_m "molar concentration at the wall of RO membrane";
    Real J_v(unit="m/s") "Solvent (water) flux in m3/m2.s";
    Real J_v1(unit="mol/(m2.s)") "Solvent (water) flux in mol/(m2.s)";
    Real J_s(unit="mol/(m2.s)") "Solute (salt) flux in mol/(m2.s)";
    Modelica.SIunits.Velocity U_b "bulk (or feed) velocity";
    Real Sc "Schimdt dimensionless number";
    Real Pe "Peclet dimensionless number L-U/D";
    Real k "mass transfer coefficient (dependent on membrane property and also function of flow conditions on the feed side)";
    Modelica.SIunits.Length D_h "hydraulic diameter, used in calculation of pressure drop";
    Real k_fb "???";
    Modelica.SIunits.Pressure delta_p "Feed pressure drop";
    Modelica.SIunits.Area A_dx "Area of partial element";
    Real error_Jv "relative error of solvent flux. This is for checking if assumption of Eq. 6 is valid";
    Modelica.SIunits.MolarDensity error_cp "error molar permeate concentration ???";
    Real osmosis(unit="m/s") "osmotic solvent (water) flux in m3/m2.s";
    Real rosmosis(unit="m/s") "reverse osmotic solvent (water) flux in m3/m2.s";

    Port_f pin_n annotation(Placement(visible = true, transformation(origin={-100,0},    extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={-100,0},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Port_f pin_n1 annotation(Placement(visible = true, transformation(origin={0,100},       extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={0,100},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Port_f pin_n4 annotation(Placement(visible = true, transformation(origin={0,-100},   extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={0,-100},   extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Port_f pin_n5 annotation(Placement(visible = true, transformation(origin={98,0},       extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={98,0},      extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    pin_n.p = p_fin;
    pin_n.Q = Q_fin;
    pin_n.c = c_fin;

    pin_n1.p = p_pin;
    pin_n1.Q = Q_pin;
    pin_n1.c = c_pin;

    pin_n4.p = p_fout;
    pin_n4.Q = -Q_fout;//-ve
    pin_n4.c = c_fout;

    pin_n5.p = p_pout;
    pin_n5.Q = -Q_pout;//-ve
    pin_n5.c = c_pout;

  //state variables correlations (density, viscosity; SI)
    rho_fin = 0.037*c_fin + 1000.3;
    rho_pin =0.037*c_pin + 1000.3;
    mu_fin = 9 * 10 ^ (-4);
  //(4.284*10^(-5)+1/(1.57*10^(-1)*(T-273.15+6.5*10^(1))^2-9.13*10^(1)))*(1+(1.54+1.998*10^(-2)*(T-273)-9.52*10^(-5)*(T-273)^2)*(c_fin*M_s/rho_fin)+(7.974-7.56*10^(-2)*(T-273)+4.724*10^(-4)*(T-273)^2)*(c_fin*M_s/rho_fin)^2);
    rho_fout = 0.037*c_fout + 1000.3;
    rho_pout = 0.037*c_pout + 1000.3;
    A_dx = dx * W;
  //mass transfer coefficient
    Sc = mu_fin / (D * rho_fin);
    U_b = (Q_fin / (2 * h_b * W)); // factor 2 for the fact that 2RO membranes speeze the permeate collecting material per leaf
  //bulk (or feed) velocity
    Pe = 2 * h_b * U_b / D;
  //Peclet dimensionless number
  //correlation from Senthilmurugan source
  k=0.55e-4;
    //k = 0.753 * (K / (2 - K)) ^ 0.5 * (D / h_b) * Sc ^ (-1 / 6) *  Modelica.Fluid.Utilities.regRoot(Pe * h_b / L_mix,0.1);
  //m/s
  // solution
  //RO equations
   J_v = A_w * (p_fin-p_pout - gamma * T * (c_m - c_pout));
      error_Jv=(A_w*(p_fin-p_pout-gamma*T*c_m)-J_v)/J_v;
       osmosis=A_w*gamma * T * (c_m - c_pout);
       rosmosis=A_w*(p_fin-p_pout);
  //m3/m2-s
    J_v1 = J_v * rho_fin / M_w;
  //mol/m2-s
    J_s = B_s * (c_m - c_pout);
  //mol/m2-s
   c_m = c_pout + (c_fin - c_pout) * exp(J_v / k);
   error_cp= c_m - (c_pout + (c_fin - c_pout) * exp(J_v / k));
  //mol/m3
  //mass balance on permeate side
    Q_pout * rho_pout = Q_pin * rho_pin + n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
  //solution balance
    Q_pout * c_pout = Q_pin * c_pin + n_leaves * 2 * A_dx * J_s;
  //salt balance
  //mass balance on feed side
    Q_fout * rho_fout = Q_fin * rho_fin -n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
  //solution balance
    Q_fout * c_fout = Q_fin * c_fin - n_leaves * 2 * A_dx * J_s;
  //salt balance
  // finding new p_f
    D_h = 2 * h_b;
  //hydrulic diameter
    //Re = rho_fin * abs(U_b) * D_h / mu_fin;
    //f = k_f * (abs(Re))^(-0.5);
    //delta_p = 2 * f * dx * U_b ^ 2 * rho_fin / D_h / 1e5;
    k_fb=18.3673e8;
    delta_p=k_fb*mu_fin*abs(U_b)^1.2;
  //bar
    if (p_fin - delta_p)>0 then
        p_fout = p_fin - delta_p;
    else
     p_fout=0;
     end if;
    p_pout = p_pin;//-delta_p*(Q_pout/Q_fout)^2;
  //Power=Q_fin*p_fin/Q_pout /(36); //kWh/m3 of desalted water
    annotation (
      conversion(noneFromVersion=""),
      Diagram(graphics={
          Polygon(
            points={{-100,100},{100,100},{100,-100},{-100,100}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-38,84},{108,0}},
            lineColor={255,255,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            fontSize=72,
            textString="Permeate"),
          Text(
            extent={{-102,4},{4,-66}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="Feed",
            fontSize=72),
          Line(
            points={{-100,100},{-100,-100},{100,-100},{100,100},{-100,100}},
            color={0,0,0},
            thickness=0.5)}),
      Icon(graphics={
          Polygon(
            points={{-100,100},{100,100},{100,-100},{-100,100}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-102,4},{4,-66}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="Feed",
            fontSize=72),
          Text(
            extent={{-38,86},{108,2}},
            lineColor={255,255,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            fontSize=72,
            textString="Permeate"),
          Line(
            points={{-100,100},{-100,-100},{100,-100},{100,100},{-100,100}},
            color={0,0,0},
            thickness=0.5)}));
  end RO_middle;

  model RO_first "First element for modeling of membrane"
    Real M_Na(unit="kg/mol") = 22.9898e-3;
    Real M_Cl(unit="kg/mol") = 35.4527e-3;
    Real M_s(unit="kg/mol") = M_Na + M_Cl;
    Real M_w(unit="kg/mol") = 18.0153e-3;
    Real D(unit="m2/s") = 1.6e-9;
    parameter Real T(unit="K") = 273 + 25;
    Real gamma(unit="Pa.m3/(K.mol)") = 8.3086;
    Real k_f = 9.6;
    parameter Real h_b(unit="m") = 0.7e-3;
    parameter Real dx(unit="m") = 0.88 / 5;
    parameter Real W( unit="m")= 1.43;
    parameter Real n_leaves = 3;
    parameter Real c_fin(unit="mol/m3")=600;
    Real L_mix(unit="m") = 0.006;
    Real K = 0.5;
    Real A_w(unit="m/(Pa.s)") = 2.08e-12;
    Real B_s(unit="m/s") = 1.11e-7;

    Real Q_fin(unit="m^3/s"), p_fin(unit="bar"),rho_fin(unit="kg/m3"), mu_fin(unit="kg/(m.s)"), p_fout(unit="bar"), c_fout(unit="mol/m3"), Q_fout(unit="m3/s"), rho_fout(unit="kg/m3"), c_pout(unit="mol/m3"), p_pout(unit="bar"), Q_pout(unit="m3/s"), rho_pout(unit="kg/m3"), c_m(unit="mol/m3"), J_v(unit="m/s"), J_v1(unit="mol/(m2.s)"), J_s(unit="mol/(m2.s)"), U_b(unit="m/s"), Sc, Pe, k, D_h(unit="m"), k_fb, delta_p(unit="bar"), A_dx(unit="m2"), error_Jv, osmosis, rosmosis;
    Port_f pin_n4 annotation(Placement(visible = true, transformation(origin={0,-100},   extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={0,-100},   extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Port_f pin_n5 annotation(Placement(visible = true, transformation(origin={100,0},      extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={100,0},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    //pin_n.p = p_fin; //declare them as eqaul to fluid ports in the RO
    //pin_n.Q = Q_fin; //declare them as eqaul to fluid ports in the RO
   // pin_n.c = c_fin;

    pin_n4.p = p_fout;
    pin_n4.Q = -Q_fout;//-ve
    pin_n4.c = c_fout;
    pin_n5.p = p_pout;
    pin_n5.Q = -Q_pout;//-ve

    pin_n5.c = c_pout;
    p_pout=0;
  //state variables correlations (density, viscosity; SI)
    rho_fin = 0.037*c_fin + 1000.3;
    mu_fin = 9 * 10 ^ (-4);
  //(4.284*10^(-5)+1/(1.57*10^(-1)*(T-273.15+6.5*10^(1))^2-9.13*10^(1)))*(1+(1.54+1.998*10^(-2)*(T-273)-9.52*10^(-5)*(T-273)^2)*(c_fin*M_s/rho_fin)+(7.974-7.56*10^(-2)*(T-273)+4.724*10^(-4)*(T-273)^2)*(c_fin*M_s/rho_fin)^2);
    rho_fout = 0.037*c_fout + 1000.3;
    rho_pout =  0.037*c_pout + 1000.3;
    A_dx = dx * W;
  //mass transfer coefficient
    Sc = mu_fin / (D * rho_fin);
    U_b = (Q_fin / (2 * h_b * W));
  //bulk (or feed) velocity
    Pe = 2 * h_b * U_b / D;
  //Peclet dimensionless number
  //correlation from Senthilmurugan source
  k=0.55e-4;
    //k = 0.753 * (K / (2 - K)) ^ 0.5 * (D / h_b) * Sc ^ (-1 / 6) * Modelica.Fluid.Utilities.regRoot(Pe * h_b / L_mix,0.1);
  //m/s
  // solution
  //RO equations
    J_v = rosmosis-osmosis;
    rosmosis=A_w*(p_fin - p_pout);
    osmosis=A_w*gamma * T * (c_m - c_pout);
  //m3/m2-s
    J_v1 = J_v * rho_fin / M_w;
  //mol/m2-s
    J_s = B_s * (c_m - c_pout);
  //mol/m2-s
    c_m = c_pout + (c_fin - c_pout) * exp(J_v / k);
    error_Jv=(A_w*(p_fin-p_pout-gamma*T*c_m)-J_v)/J_v;
  //mol/m3
  //mass balance on permeate side
    Q_pout * rho_pout = n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
  //solution balance
    Q_pout * c_pout = n_leaves * 2 * A_dx * J_s;
  //salt balance
  //mass balance on feed side
    Q_fout * rho_fout = Q_fin * rho_fin -n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
  //solution balance
    Q_fout * c_fout = Q_fin * c_fin -n_leaves * 2 * A_dx * J_s;
  //salt balance
  // finding new p_f
    D_h = 2 * h_b;
  //hydrulic diameter
    //Re = rho_fin * abs(U_b) * D_h / mu_fin;
    //f = k_f * (Re)^(-0.5);
    //delta_p = 2 * f * dx * U_b ^ 2 * rho_fin / D_h / 1e5;
    k_fb=18.3673e8;
    delta_p=k_fb*mu_fin*abs(U_b)^1.2;
  //bar
  if (p_fin - delta_p)>0 then
        p_fout = p_fin - delta_p;
    else
     p_fout=0;
     end if;
  //Power=Q_fin*p_fin/Q_pout /(36); //kWh/m3 of desalted water
    annotation (                                Diagram(graphics={
          Polygon(
            points={{-100,100},{100,100},{100,-100},{-100,100}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-38,84},{108,0}},
            lineColor={255,255,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            fontSize=72,
            textString="Permeate"),
          Text(
            extent={{-102,4},{4,-66}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="Feed",
            fontSize=72),
          Line(
            points={{-100,100},{-100,-100},{100,-100},{100,100},{-100,100}},
            color={0,0,0},
            thickness=0.5)}),
      Icon(graphics={
          Polygon(
            points={{-100,100},{100,100},{100,-100},{-100,100}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-38,86},{108,2}},
            lineColor={255,255,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            fontSize=72,
            textString="Permeate"),
          Text(
            extent={{-102,4},{4,-66}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="Feed",
            fontSize=72),
          Line(
            points={{-100,100},{-100,-100},{100,-100},{100,100},{-100,100}},
            color={0,0,0},
            thickness=0.5)}));
  end RO_first;

  model RO_last "Last element for modeling of membrane"
    Real M_Na(unit="kg/mol") = 22.9898e-3;
    Real M_Cl(unit="kg/mol") = 35.4527e-3;
    Real M_s(unit="kg/mol") = M_Na + M_Cl;
    Real M_w(unit="kg/mol") = 18.0153e-3;
    Real D(unit="m2/s") = 1.6e-9;
    parameter Real T(unit="K") = 273 + 25;
    Real gamma(unit="Pa.m3/(K.mol)") = 8.3086;
    Real k_f = 9.6;
    parameter Real h_b(unit="m") = 0.7e-3;
    parameter Real dx(unit="m") = 0.88 / 5;
    parameter Real W( unit="m")= 1.43;
    parameter Real n_leaves = 3;
    Real L_mix(unit="m") = 0.006;
    Real K = 0.5;
    Real A_w(unit="m/(Pa.s)") = 2.08e-12;
    Real B_s(unit="m/s") = 1.11e-7;

  Real c_fin(unit="mol/m3"), p_fin(unit="bar"), Q_fin(unit="m3/s"), rho_fin(unit="kg/m3"), mu_fin(unit="kg/(m.s)"), p_fout(unit="bar"), c_fout(unit="mol/m3"), Q_fout(unit="m3/s"), rho_fout(unit="kg/m3"), c_pin(unit="mol/m3"), p_pin(unit="bar"), Q_pin(unit="m3/s"), rho_pin(unit="kg/m3"), c_pout(unit="mol/m3"), p_pout(unit="bar"), Q_pout(unit="m3/s"), rho_pout(unit="kg/m3"), c_m(unit="mol/m3"), J_v(unit="m/s"), J_v1(unit="mol/(m2.s)"), J_s(unit="mol/(m2.s)"), U_b(unit="m/s"), Sc, Pe, k, D_h(unit="m"), k_fb, delta_p(unit="bar"), A_dx(unit="m2"), error_Jv, error_cp, osmosis, rosmosis;

    Port_f pin_n annotation(Placement(visible = true, transformation(origin={-100,0},    extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={-100,0},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Port_f pin_n1 annotation(Placement(visible = true, transformation(origin={0,98},        extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin={0,100},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    pin_n.p = p_fin;
    pin_n.Q = Q_fin;
    pin_n.c = c_fin;

    pin_n1.p = p_pin;
    pin_n1.Q = Q_pin;
    pin_n1.c = c_pin;

    //pin_n2.Q = Q_fout;
    //pin_n2.c = c_fout;

  //state variables correlations (density, viscosity; SI)
    rho_fin = 0.037 * c_fin + 1000.3;
    rho_pin = 0.037 * c_pin + 1000.3;
    mu_fin = 9 * 10 ^ (-4);
  //(4.284*10^(-5)+1/(1.57*10^(-1)*(T-273.15+6.5*10^(1))^2-9.13*10^(1)))*(1+(1.54+1.998*10^(-2)*(T-273)-9.52*10^(-5)*(T-273)^2)*(c_fin*M_s/rho_fin)+(7.974-7.56*10^(-2)*(T-273)+4.724*10^(-4)*(T-273)^2)*(c_fin*M_s/rho_fin)^2);
    rho_fout = 0.037 * c_fout + 1000.3;
    rho_pout = 0.037 * c_pout + 1000.3;
    A_dx = dx * W;
  //mass transfer coefficient
    Sc = mu_fin / (D * rho_fin);
    U_b = Q_fin / (2 * h_b * W);
  //bulk (or feed) velocity
    Pe = 2 * h_b * U_b / D;
  //Peclet dimensionless number
  //correlation from Senthilmurugan source
  k=0.55e-4;
    //k = 0.753 * (K / (2 - K)) ^ 0.5 * (D / h_b) * Sc ^ (-1 / 6) *  Modelica.Fluid.Utilities.regRoot(Pe * h_b / L_mix,0.1);
  //m/s
  // solution
  //RO equations
    J_v = A_w * (p_fin - p_pout - gamma * T * (c_m - c_pout));
    error_Jv = (A_w * (p_fin - p_pout - gamma * T * c_m) - J_v) / J_v;
     osmosis=A_w*gamma * T * (c_m - c_pout);
     rosmosis=A_w * (p_fin - p_pout);
  //m3/m2-s
    J_v1 = J_v * rho_fin / M_w;
  //mol/m2-s
    J_s = B_s * (c_m - c_pout);
  //mol/m2-s
    c_m = c_pout + (c_fin - c_pout) * exp(J_v / k);
     error_cp= c_m - (c_pout + (c_fin - c_pout) * exp(J_v / k));
  //mol/m3
  //mass balance on permeate side
    Q_pout * rho_pout = Q_pin * rho_pin + n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
  //solution balance
    Q_pout * c_pout = Q_pin * c_pin + n_leaves * 2 * A_dx * J_s;
  //salt balance
  //mass balance on feed side
    Q_fout * rho_fout = Q_fin * rho_fin - n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
  //solution balance
    Q_fout * c_fout = Q_fin * c_fin - n_leaves * 2 * A_dx * J_s;
  //salt balance
  // finding new p_f
    D_h = 2 * h_b;
  //hydrulic diameter
    //Re = rho_fin * abs(U_b) * D_h / mu_fin;
    //f = k_f * (abs(Re))^(-0.5);
    //delta_p = 2 * f * dx * U_b ^ 2 * rho_fin / D_h / 1e5;
    k_fb=18.3673e8;
    delta_p=k_fb*mu_fin*abs(U_b)^1.2;
  //bar
    if (p_fin - delta_p)>0 then
        p_fout = p_fin - delta_p;
    else
     p_fout=0;
     end if;
    p_pout = p_pin;
  //-delta_p*(Q_pout/Q_fout)^2;
  //Power=Q_fin*p_fin/Q_pout /(36); //kWh/m3 of desalted water
    annotation (
      Diagram(graphics={
          Polygon(
            points={{-100,100},{100,100},{100,-100},{-100,100}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-38,84},{108,0}},
            lineColor={255,255,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            fontSize=72,
            textString="Permeate"),
          Text(
            extent={{-102,4},{4,-66}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="Feed",
            fontSize=72),
          Line(
            points={{-100,100},{-100,-100},{100,-100},{100,100},{-100,100}},
            color={0,0,0},
            thickness=0.5)}),
      Icon(graphics={
          Polygon(
            points={{-100,100},{100,100},{100,-100},{-100,100}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-38,86},{108,2}},
            lineColor={255,255,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            fontSize=72,
            textString="Permeate"),
          Text(
            extent={{-102,4},{4,-66}},
            lineColor={0,0,255},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="Feed",
            fontSize=72),
          Line(
            points={{-100,100},{-100,-100},{100,-100},{100,100},{-100,100}},
            color={0,0,0},
            thickness=0.5)}));
  end RO_last;

  connector Port_f "fluid pin for RO"
    Real p(unit="bar"),c(unit="mol/m3");
     flow Real Q(unit="m3/s");
    annotation(defaultComponentName = "pin_n", Documentation(info = "<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>",   revisions = "<html>
<dl>
<dt><i>1998</i></dt>
<dd>by Christoph Clauss initially implemented
</dd>
</dl>
</html>"),   Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid), Text(extent={{
                -96,100},{104,40}},                                                                                                                                                                                                        lineColor=
                {0,0,255},
            textString="RO port")}));
  end Port_f;

  model RO_array_test

    parameter Integer n=2 "number of modules";
    RO_middle module[n]
    annotation (Placement(transformation(
    extent={{-10,-10},{10,10}},
    origin={-30,20})));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end RO_array_test;
  annotation (uses(
      OpenHydraulics(version="1.0"),
      Modelica(version="3.2.2"),
      RO_middle(version="1")));
end RO_hydraulics;
