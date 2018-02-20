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
    connect(cosine1.y, source.control) annotation(Line(points={{-107,8},{-68,8},
            {-68,8},{-68,8}},                                                                              color = {0, 0, 127}));
    connect(source2.port, rO1.port_fout)
      annotation (Line(points={{-22,-18},{-22,8.2}}, color={255,0,0}));
    connect(source.port, rO1.port_fin)
      annotation (Line(points={{-58,18},{-32,18}}, color={255,0,0}));
  end RO_system;

  model RO "Series of segments of one RO module modeled as one element"

  import OpenHydraulics.Interfaces;
    extends OpenHydraulics.Interfaces.PartialFluidComponent;
    Modelica.SIunits.EnergyDensity SpecificEnergy;
    RO_middle rO_middle1(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin={-40,10},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle rO_middle2(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin={0,10},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle rO_middle3(dx = 0.88 / 5) annotation(Placement(visible = true, transformation(origin={38,10},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    OpenHydraulics.Interfaces.FluidPort port_fin "feed inlet port" annotation (
        Placement(
        visible=true,
        transformation(
          origin={-100,0},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={-100,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    OpenHydraulics.Interfaces.FluidPort port_fout "feed outlet port (reject)"
      annotation (Placement(
        visible=true,
        transformation(
          origin={100,0},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={0,-98},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    RO_middle rO_first(dx=0.88/5) annotation (Placement(visible=true,
          transformation(
          origin={-74,10},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    parameter Modelica.SIunits.MolarDensity c_s=600 "salt water inlet molar salt concentration";

    RO_middle rO_last(dx=0.88/5) annotation (Placement(visible=true,
          transformation(
          origin={74,10},
          extent={{-10,-10},{10,10}},
          rotation=0)));
  equation
    connect(rO_middle2.port_pout, rO_middle3.port_pin) annotation (Line(points={{9.8,10},
            {18,10},{18,28},{38,28},{38,20}},             color={0,0,255}));
    connect(rO_middle2.port_fout, rO_middle3.port_fin) annotation (Line(points={{0,0},{0,
            -10},{22,-10},{22,10},{28,10}},               color={0,0,255}));
    connect(rO_middle1.port_pout, rO_middle2.port_pin) annotation (Line(points={{-30.2,
            10},{-20,10},{-20,28},{0,28},{0,20}},   color={0,0,255}));
    connect(rO_middle1.port_fout, rO_middle2.port_fin) annotation (Line(points={{-40,0},
            {-40,-10},{-16,-10},{-16,10},{-10,10}},color={0,0,255}));
    // set constratints to first element
    rO_first.port_fin.c = c_s;
    rO_first.port_fin.p = port_fin.p;
    rO_first.port_fin.Q = port_fin.m_flow/oil.density(port_fin.p);

    // set constraints to last element
    rO_last.port_fout.p = port_fout.p;
    -rO_last.port_fout.Q = port_fout.m_flow/oil.density(port_fout.p);
    rO_last.port_pout.p = 1.025E5;

    SpecificEnergy=rO_first.port_fin.Q*rO_first.port_fout.p/(rO_last.port_pout.Q)/(3.6e6);

    connect(rO_first.port_pout, rO_middle1.port_pin) annotation (Line(points={{-64.2,
            10},{-58,10},{-58,28},{-40,28},{-40,20}},          color={0,0,255}));
    connect(rO_middle1.port_fin, rO_first.port_fout) annotation (Line(points={{-50,10},
            {-54,10},{-54,-10},{-74,-10},{-74,0}},     color={0,0,255}));
    connect(rO_middle3.port_fout, rO_last.port_fin) annotation (Line(points={{38,0},
            {38,-10},{58,-10},{58,10},{64,10}}, color={0,0,255}));
    connect(rO_middle3.port_pout, rO_last.port_pin) annotation (Line(points={{47.8,
            10},{52,10},{52,30},{74,30},{74,20}}, color={0,0,255}));
    annotation (                                                                Icon(graphics={  Polygon(lineColor = {0, 0, 255}, fillColor = {7, 0, 255},
              fillPattern =                                                                                                                                              FillPattern.Solid, points = {{-100, 100}, {100, 100}, {100, -100}, {100, -100}, {-100, 100}}), Line(points = {{-100, -100}, {-100, 100}, {100, 100}, {100, -100}, {-100, -100}})}));
  end RO;

  model RO_middle "middle element of RO membrane"
  protected
    type MolarFlux = Real (final quantity="MolarFlux", final unit="mol/(m2.s)");
    type VolumetricFlux = Real (final quantity="VolumetricFlux", final unit="m3/(m2.s)");
    parameter Modelica.SIunits.MolarMass M_Na = 22.9898e-3 "Molar mass of sodium [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_Cl = 35.4527e-3 "Molar mass of chlorine [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_s = M_Na + M_Cl "Molar mass of sodium chloride [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_w = 18.0153e-3 "Molar mass of water [kg/mol]";
    parameter Modelica.SIunits.DiffusionCoefficient D = 1.6e-9 "Water - Salt mass diffusivity https://en.wikipedia.org/wiki/Mass_diffusivity [m2/s]";
    parameter Modelica.SIunits.Temperature T = 273 + 25;
    parameter Modelica.SIunits.MolarHeatCapacity gamma = 8.3086 "Molar gas constant R [J/mol.K]";
    parameter Modelica.SIunits.DynamicViscosity mu_fin = 9e-4 "feed inlet dynamic viscosity";


    // RO module specific parameters
  public
    parameter Real k_f = 9.6 "friction factor, find experimentally";
    parameter Modelica.SIunits.Height h_b = 0.7e-3 "height of feed channel spacers";
    parameter Modelica.SIunits.Length dx = 0.88 / 5 "Length of finite element of RO module";
    parameter Modelica.SIunits.Length W = 1.43 "width of element";
    parameter Modelica.SIunits.Area A_dx = dx * W "Area of finite element of RO module";
    parameter Integer n_leaves = 3 "number of leaves";
    parameter Modelica.SIunits.Length L_mix = 0.006 "characteristic length of mixing net (spacer)";
    parameter Real K = 0.5 "mixing coefficient of net";
    parameter Real A_w(unit="m/(Pa.s)") = 2.08e-12 "RO membrane solvent permeability coefficient";
    parameter Real B_s(unit="m/s") = 1.11e-7 "membrane coefficient to describe salt diffusion as a function of concentration difference, membrane specific";
  protected
    parameter Modelica.SIunits.Length D_h = 2 * h_b "Hydraulic diameter of gap on feed side";

  public
    Modelica.SIunits.MolarDensity c_fin "feed inlet molar salt concentration";
    Modelica.SIunits.Pressure p_fin "feed inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_fin "feed inlet volume flow rate";
    Modelica.SIunits.Density rho_fin "feed inlet density";

    Modelica.SIunits.Pressure p_fout "feed outlet pressure";
    Modelica.SIunits.MolarDensity c_fout "feed outlet molar salt concentration";
    Modelica.SIunits.VolumeFlowRate Q_fout "feed outlet volume flow rate";
    Modelica.SIunits.Density rho_fout "feed outlet density";

    Modelica.SIunits.MolarDensity c_pin "permeate inlet molar salt concentration";
    Modelica.SIunits.Pressure p_pin "permeate inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pin "permeate inlet volume flow rate";
    Modelica.SIunits.Density rho_pin "permeate inlet density";

    Modelica.SIunits.MolarDensity c_pout "pearmeate outlet molar salt concentration";
    Modelica.SIunits.Pressure p_pout "permeate outlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pout "permeate outlet volume flow rate";
    Modelica.SIunits.Density rho_pout "permeate outlet density";

    Modelica.SIunits.MolarDensity c_m "molar salt concentration on wall of RO membrane";
    VolumetricFlux J_v "Molar solvent (water) flux through membrane in [m3/(m2.s)]";
    MolarFlux J_v1 "Volumetric solvent (water) flux through membrane in [mol/(m2.s)]";
    MolarFlux J_s "Molar salt flux through membrane [mol/(m2.s)]";
    Modelica.SIunits.Velocity U_b "average feed velocity (based on inlet flow, spacer height and width of membrane)";
    Real Sc "Dimensionless Schmidt number";
    Real Pe "Dimensionless Peclet number";
    Real k "Mass transfer coefficient, dependent on membrane properties and flow conditions of feed side";
    Real k_fb "Loss factor for laminar turbulent flow in feed channel";
    Modelica.SIunits.Pressure delta_p "Pressure loss in feed channel";
    Real error_Jv "Relative error of solvent flux to check if applied assumptions are fine";
    Modelica.SIunits.MolarDensity error_cp "error of molar permeate salt concentration";
    VolumetricFlux osmosis "osmotic flux. Given the concentration difference at feed wall and permeate, this is how much is going through the membrane in osmotic direction";
    VolumetricFlux rosmosis "reverse osmotic flux. Given the pressure on feed and permeate sides, this is how much is going throught the membrane in reverse osmotic direction";

    Port_f port_fin "feed inlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={-100,0},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={-100,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Port_f port_pin "permeate inlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={0,100},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={0,100},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Port_f port_fout "feed outlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Port_f port_pout "permeate outlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={98,0},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={98,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));
  equation
    port_fin.p = p_fin;
    port_fin.Q = Q_fin;
    port_fin.c = c_fin;

    port_pin.p = p_pin;
    port_pin.Q = Q_pin;
    port_pin.c = c_pin;

    port_fout.p = p_fout;
    port_fout.Q = -Q_fout; //-ve
    port_fout.c = c_fout;

    port_pout.p = p_pout;
    port_pout.Q = -Q_pout; //-ve
    port_pout.c = c_pout;

  //state variables correlations (density, viscosity; SI)
    rho_fin = 0.037*c_fin + 1000.3;
    rho_pin = 0.037*c_pin + 1000.3;

  //(4.284*10^(-5)+1/(1.57*10^(-1)*(T-273.15+6.5*10^(1))^2-9.13*10^(1)))*(1+(1.54+1.998*10^(-2)*(T-273)-9.52*10^(-5)*(T-273)^2)*(c_fin*M_s/rho_fin)+(7.974-7.56*10^(-2)*(T-273)+4.724*10^(-4)*(T-273)^2)*(c_fin*M_s/rho_fin)^2);
    rho_fout = 0.037*c_fout + 1000.3;
    rho_pout = 0.037*c_pout + 1000.3;

  //mass transfer coefficient
    Sc = mu_fin / (D * rho_fin);
    U_b = (Q_fin / (2 * h_b * W));  //bulk (or feed) velocity
    Pe = 2 * h_b * U_b / D;         //Peclet dimensionless number

  //correlation from Senthilmurugan source
    k=0.55e-4;
    //k = 0.753 * (K / (2 - K)) ^ 0.5 * (D / h_b) * Sc ^ (-1 / 6) *  Modelica.Fluid.Utilities.regRoot(Pe * h_b / L_mix,0.1);
  //m/s

  // RO equations
   J_v = A_w * (p_fin-p_pout - gamma * T * (c_m - c_pout));
       error_Jv=(A_w*(p_fin-p_pout-gamma*T*c_m)-J_v)/J_v;
       osmosis=A_w*gamma * T * (c_m - c_pout);
       rosmosis=A_w*(p_fin-p_pout);
   J_v1 = J_v * rho_fin / M_w;
   J_s = B_s * (c_m - c_pout);
   c_m = c_pout + (c_fin - c_pout) * exp(J_v / k);
   error_cp= c_m - (c_pout + (c_fin - c_pout) * exp(J_v / k));

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
    //Re = rho_fin * abs(U_b) * D_h / mu_fin;
    //f = k_f * (abs(Re))^(-0.5);
    //delta_p = 2 * f * dx * U_b ^ 2 * rho_fin / D_h / 1e5;
    k_fb=18.3673e8;
    delta_p=k_fb*mu_fin*abs(U_b)^1.2;

    if (p_fin - delta_p)>0 then
        p_fout = p_fin - delta_p;
    else
      p_fout=0;
    end if;
  // deltap on permeate side
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
    Modelica.SIunits.Pressure p "pressure";
    Modelica.SIunits.MolarDensity c "Molar salt concentration";
    flow Modelica.SIunits.VolumeFlowRate Q "volumetric flow rate";
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

  model RO_middle_stream
    "middle element of RO membrane with stream connectors"
  protected
    type MolarFlux = Real (final quantity="MolarFlux", final unit="mol/(m2.s)");
    type VolumetricFlux = Real (final quantity="VolumetricFlux", final unit="m3/(m2.s)");
    parameter Modelica.SIunits.MolarMass M_Na = 22.9898e-3 "Molar mass of sodium [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_Cl = 35.4527e-3 "Molar mass of chlorine [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_s = M_Na + M_Cl "Molar mass of sodium chloride [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_w = 18.0153e-3 "Molar mass of water [kg/mol]";
    parameter Modelica.SIunits.DiffusionCoefficient D = 1.6e-9 "Water - Salt mass diffusivity https://en.wikipedia.org/wiki/Mass_diffusivity [m2/s]";
    parameter Modelica.SIunits.Temperature T = 273 + 25;
    parameter Modelica.SIunits.MolarHeatCapacity gamma = 8.3086 "Molar gas constant R [J/mol.K]";
    // RO module specific parameters
  public
    parameter Real k_f = 9.6 "friction factor, find experimentally";
    parameter Modelica.SIunits.Height h_b = 0.7e-3 "height of feed channel spacers";
    parameter Modelica.SIunits.Length dx = 0.88 / 5 "Length of finite element of RO module";
    parameter Modelica.SIunits.Length W = 1.43 "width of element";
    parameter Integer n_leaves = 3 "number of leaves";
    parameter Modelica.SIunits.Length L_mix = 0.006 "characteristic length of mixing net (spacer)";
    parameter Real K = 0.5 "mixing coefficient of net";
    parameter Real A_w(unit="m/(Pa.s)") = 2.08e-12 "RO membrane solvent permeability coefficient";
    Real B_s(unit="m/s") = 1.11e-7 "membrane coefficient to describe salt diffusion as a function of concentration difference, membrane specific";

    Modelica.SIunits.MolarDensity c_fin "feed inlet molar salt concentration";
    Modelica.SIunits.Pressure p_fin "feed inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_fin "feed inlet volume flow rate";
    Modelica.SIunits.Density rho_fin "feed inlet density";
    Modelica.SIunits.DynamicViscosity mu_fin "feed inlet dynamic viscosity";

    Modelica.SIunits.Pressure p_fout "feed outlet pressure";
    Modelica.SIunits.MolarDensity c_fout "feed outlet molar salt concentration";
    Modelica.SIunits.VolumeFlowRate Q_fout "feed outlet volume flow rate";
    Modelica.SIunits.Density rho_fout "feed outlet density";

    Modelica.SIunits.MolarDensity c_pin "permeate inlet molar salt concentration";
    Modelica.SIunits.Pressure p_pin "permeate inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pin "permeate inlet volume flow rate";
    Modelica.SIunits.Density rho_pin "permeate inlet density";

    Modelica.SIunits.MolarDensity c_pout "pearmeate outlet molar salt concentration";
    Modelica.SIunits.Pressure p_pout "permeate outlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pout "permeate outlet volume flow rate";
    Modelica.SIunits.Density rho_pout "permeate outlet density";

    Modelica.SIunits.MolarDensity c_m "molar salt concentration on wall of RO membrane";
    VolumetricFlux J_v "Molar solvent (water) flux through membrane in [m3/(m2.s)]";
    MolarFlux J_v1 "Volumetric solvent (water) flux through membrane in [mol/(m2.s)]";
    MolarFlux J_s "Molar salt flux through membrane [mol/(m2.s)]";
    Modelica.SIunits.Velocity U_b "average feed velocity (based on inlet flow, spacer height and width of membrane)";
    Real Sc "Dimensionless Schmidt number";
    Real Pe "Dimensionless Peclet number";
    Real k "Mass transfer coefficient, dependent on membrane properties and flow conditions of feed side";
    Modelica.SIunits.Length D_h "Hydraulic diameter of gap on feed side";
    Real k_fb "Loss factor for laminar turbulent flow in feed channel";
    Modelica.SIunits.Pressure delta_p "Pressure loss in feed channel";
    Modelica.SIunits.Area A_dx "Area of finite element of RO module";
    Real error_Jv "Relative error of solvent flux to check if applied assumptions are fine";
    Modelica.SIunits.MolarDensity error_cp "error of molar permeate salt concentration";
    VolumetricFlux osmosis "osmotic flux. Given the concentration difference at feed wall and permeate, this is how much is going through the membrane in osmotic direction";
    VolumetricFlux rosmosis "reverse osmotic flux. Given the pressure on feed and permeate sides, this is how much is going throught the membrane in reverse osmotic direction";

    Port_f port_fin "feed inlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={-100,0},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={-100,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Port_f port_pin "permeate inlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={0,100},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={0,100},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Port_f port_fout "feed outlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Port_f port_pout "permeate outlet port" annotation (Placement(
        visible=true,
        transformation(
          origin={98,0},
          extent={{-10,-10},{10,10}},
          rotation=0),
        iconTransformation(
          origin={98,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));
  equation
    port_fin.p = p_fin;
    port_fin.Q = Q_fin;
    port_fin.c = c_fin;

    port_pin.p = p_pin;
    port_pin.Q = Q_pin;
    port_pin.c = c_pin;

    port_fout.p = p_fout;
    port_fout.Q = -Q_fout; //-ve
    port_fout.c = c_fout;

    port_pout.p = p_pout;
    port_pout.Q = -Q_pout; //-ve
    port_pout.c = c_pout;

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
    U_b = (Q_fin / (2 * h_b * W));
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
  end RO_middle_stream;
  annotation (uses(
      OpenHydraulics(version="1.0"),
      Modelica(version="3.2.2"),
      RO_middle(version="1")));
end RO_hydraulics;
