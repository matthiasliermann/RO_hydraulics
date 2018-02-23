package RO_hydraulics "package for modeling reverse osmosis membrane process"
  model RO_middle_stream "middle element of RO membrane with stream connectors"
    extends Interfaces.PartialFourPort;
  protected
    type MolarFlux = Real(final quantity = "MolarFlux", final unit = "mol/(m2.s)");
    type VolumetricFlux = Real(final quantity = "VolumetricFlux", final unit = "m3/(m2.s)");
    parameter Modelica.SIunits.MolarMass M_Na = 22.9898e-3 "Molar mass of sodium [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_Cl = 35.4527e-3 "Molar mass of chlorine [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_s = M_Na + M_Cl "Molar mass of sodium chloride [kg/mol]";
    parameter Modelica.SIunits.MolarMass M_w = 18.0153e-3 "Molar mass of water [kg/mol]";
    parameter Modelica.SIunits.DiffusionCoefficient D = 1.6e-9 "Water - Salt mass diffusivity https://en.wikipedia.org/wiki/Mass_diffusivity [m2/s]";
    parameter Modelica.SIunits.Temperature T = 273 + 25;
    parameter Modelica.SIunits.MolarHeatCapacity gamma = 8.3086 "Molar gas constant R [J/mol.K]";
    // RO module specific parameters
    parameter Modelica.SIunits.Area A_dx = dx * W "Area of finite element of RO module";
    parameter Modelica.SIunits.DynamicViscosity mu_fin = 9e-4 "feed inlet dynamic viscosity";
    parameter Real k = 0.55e-4 "Mass transfer coefficient, dependent on membrane properties and flow conditions of feed side";
    // correlation from Senthilmurugan source
    // alternative make dependent of flow conditions, but is numerically difficult
    // k = 0.753 * (K / (2 - K)) ^ 0.5 * (D / h_b) * Sc ^ (-1 / 6) *  Modelica.Fluid.Utilities.regRoot(Pe * h_b / L_mix,0.1);
    // m/s
    parameter Modelica.SIunits.Length D_h = 2 * h_b "Hydraulic diameter of gap on feed side";
    parameter Real k_fb = 18.3673e8 "Loss factor for laminar turbulent flow in feed channel";
  public
    parameter Real k_f = 9.6 "friction factor, find experimentally";
    parameter Modelica.SIunits.Height h_b = 0.7e-3 "height of feed channel spacers";
    parameter Modelica.SIunits.Length dx = 0.88 / 5 "Length of finite element of RO module";
    parameter Modelica.SIunits.Length W = 1.43 "width of element";
    parameter Integer n_leaves = 3 "number of leaves";
    parameter Modelica.SIunits.Length L_mix = 0.006 "characteristic length of mixing net (spacer)";
    parameter Real K = 0.5 "mixing coefficient of net";
    parameter Real A_w(unit = "m/(Pa.s)") = 2.08e-12 "RO membrane solvent permeability coefficient";
    parameter Modelica.SIunits.MolarDensity c_start = 600 "start value for molar salt concentration of feed inlet";
    parameter Modelica.SIunits.Pressure p_start = 50e5 "start value for pressure at feed inlet";
  public
    Real B_s(unit = "m/s") = 1.11e-7 "membrane coefficient to describe salt diffusion as a function of concentration difference, membrane specific";
    Modelica.SIunits.MolarDensity c_fin(start = c_start) "feed inlet molar salt concentration";
    Modelica.SIunits.Pressure p_fin(start = p_start) "feed inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_fin "feed inlet volume flow rate";
    Modelica.SIunits.Density rho_fin "feed inlet density";
    Modelica.SIunits.Pressure p_fout "feed outlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_fout "feed outlet volume flow rate";
    Modelica.SIunits.Density rho_fout "feed outlet density";
    Modelica.SIunits.MolarDensity c_pin "permeate inlet molar salt concentration";
    Modelica.SIunits.Pressure p_pin "permeate inlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pin "permeate inlet volume flow rate";
    Modelica.SIunits.Density rho_pin "permeate inlet density";
    Modelica.SIunits.Pressure p_pout "permeate outlet pressure";
    Modelica.SIunits.VolumeFlowRate Q_pout "permeate outlet volume flow rate";
    Modelica.SIunits.Density rho_pout "permeate outlet density";
    VolumetricFlux J_v "Molar solvent (water) flux through membrane in [m3/(m2.s)]";
    MolarFlux J_v1 "Volumetric solvent (water) flux through membrane in [mol/(m2.s)]";
    MolarFlux J_s "Molar salt flux through membrane [mol/(m2.s)]";
    Modelica.SIunits.Velocity U_b "average feed velocity (based on inlet flow, spacer height and width of membrane)";
    Real Sc "Dimensionless Schmidt number";
    Real Pe "Dimensionless Peclet number";
    Real error_Jv "Relative error of solvent flux to check if applied assumptions are fine";
    Modelica.SIunits.MolarDensity error_cp "error of molar permeate salt concentration";
    VolumetricFlux osmosis "osmotic flux. Given the concentration difference at feed wall and permeate, this is how much is going through the membrane in osmotic direction";
    VolumetricFlux rosmosis "reverse osmotic flux. Given the pressure on feed and permeate sides, this is how much is going throught the membrane in reverse osmotic direction";
    Modelica.SIunits.Pressure dp_membrane "pressure difference across membrane";
    Modelica.SIunits.Pressure dp_osmotic "osmotic pressure across membrane";
    Modelica.SIunits.MolarDensity cp "Molar concentration of salt in permeate outlet";
    Modelica.SIunits.MolarDensity cf "Molar concentration of salt feed outlet";
    Modelica.SIunits.MolarDensity cm "Molar concentration of salt on membrane wall of feed side due to concentration polarization";
    Modelica.SIunits.Density cp_(displayUnit = "g/L") "Salt density in permeate";
    Modelica.SIunits.Density cf_(displayUnit = "g/L") "Salt density in reject";
  equation
    port_feed_a.p = p_fin;
    port_feed_a.m_flow = Q_fin * rho_fin;
    port_feed_a.c = if port_feed_a.m_flow > 0 then inStream(port_feed_a.c) else cf;
    port_feed_b.p = p_fout;
    port_feed_b.m_flow = -Q_fout * rho_fout;
    port_feed_b.c = if port_feed_b.m_flow > 0 then inStream(port_feed_b.c) else cf;
    c_fin = if port_feed_a.m_flow > 0 then port_feed_a.c else port_feed_b.c;
    dp_membrane = p_fin - p_pout;
    port_permeate_a.p = p_pin;
    port_permeate_a.m_flow = Q_pin * rho_pin;
    port_permeate_a.c = if port_permeate_a.m_flow > 0 then inStream(port_permeate_a.c) else cp;
    port_permeate_b.p = p_pout;
    port_permeate_b.m_flow = -Q_pout * rho_pout;
    port_permeate_b.c = if port_permeate_b.m_flow > 0 then inStream(port_permeate_b.c) else cp;
    c_pin = if port_permeate_a.m_flow > 0 then port_permeate_a.c else port_permeate_b.c;
    rho_fin = 0.037 * c_fin + 1000.3;
    //state variables correlations (density, viscosity; SI)
    rho_pin = 0.037 * c_pin + 1000.3;
    //(4.284*10^(-5)+1/(1.57*10^(-1)*(T-273.15+6.5*10^(1))^2-9.13*10^(1)))*(1+(1.54+1.998*10^(-2)*(T-273)-9.52*10^(-5)*(T-273)^2)*(c_fin*M_s/rho_fin)+(7.974-7.56*10^(-2)*(T-273)+4.724*10^(-4)*(T-273)^2)*(c_fin*M_s/rho_fin)^2);
    rho_fout = 0.037 * cf + 1000.3;
    rho_pout = 0.037 * cp + 1000.3;
    //mass transfer coefficient
    Sc = mu_fin / (D * rho_fin);
    U_b = Q_fin / (2 * h_b * W);
    //bulk (or feed) velocity
    Pe = 2 * h_b * U_b / D;
    //Peclet dimensionless number
    // solution
    //RO equations
    J_v = A_w * (dp_membrane - gamma * T * (cm - cp));
    error_Jv = (A_w * (dp_membrane - gamma * T * cm) - J_v) / J_v;
    dp_osmotic = gamma * T * (cm - cp);
    osmosis = A_w * dp_osmotic;
    rosmosis = A_w * dp_membrane;
    //m3/m2-s
    J_v1 = J_v * rho_fin / M_w;
    //mol/m2-s
    J_s = B_s * (cm - cp);
    //mol/m2-s
    // flow across membrane
    m_flow_membrane = n_leaves * 2 * A_dx * (J_v1 * M_w + J_s * M_s);
    cm = cp + (c_fin - cp) * exp(J_v / k);
    cp_ = cp*M_s;
    cf_ = cf*M_w;
    error_cp = cm - (cp + (c_fin - cp) * exp(J_v / k));
    //mol/m3
    //salt balance on permeate side
    Q_pout * cp = Q_pin * c_pin + n_leaves * 2 * A_dx * J_s;
    //salt balance on feed side
    Q_fout * cf = Q_fin * c_fin - n_leaves * 2 * A_dx * J_s;
    // finding new p_f
    // Re = rho_fin * abs(U_b) * D_h / mu_fin;
    // f = k_f * (abs(Re))^(-0.5);
    // dp_feed = 2 * f * dx * U_b ^ 2 * rho_fin / D_h / 1e5;
    dp_feed = k_fb * mu_fin * sign(U_b) * abs(U_b) ^ 1.2;
    //bar
    dp_permeate = 0;
    //-dp_feed*(Q_pout/Q_fout)^2;
    //Power=Q_fin*p_fin/Q_pout /(36); //kWh/m3 of desalted water
    annotation (
      conversion(noneFromVersion = ""),
      Diagram(graphics={  Polygon(points = {{-100, 100}, {100, 100}, {100, -100}, {-100, 100}}, lineColor = {0, 0, 255},
              lineThickness =                                                                                                            0.5, fillColor = {0, 0, 255},
              fillPattern =                                                                                                                                                          FillPattern.Solid), Text(extent = {{-38, 84}, {108, 0}}, lineColor = {255, 255, 255},
              lineThickness =                                                                                                                                                                                                        0.5, fillColor = {0, 0, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, fontSize = 72, textString = "Permeate"), Text(extent = {{-102, 4}, {4, -66}}, lineColor = {0, 0, 255},
              lineThickness =                                                                                                                                                                                                        0.5, fillColor = {0, 0, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, textString = "Feed", fontSize = 72), Line(points = {{-100, 100}, {-100, -100}, {100, -100}, {100, 100}, {-100, 100}}, color = {0, 0, 0}, thickness = 0.5)}),
      Icon(graphics={  Polygon(points = {{-100, 100}, {100, 100}, {100, -100}, {-100, 100}}, lineColor = {0, 0, 255},
              lineThickness =                                                                                                         0.5, fillColor = {0, 0, 255},
              fillPattern =                                                                                                                                                       FillPattern.Solid), Text(extent = {{-102, 4}, {4, -66}}, lineColor = {0, 0, 255},
              lineThickness =                                                                                                                                                                                                        0.5, fillColor = {0, 0, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, textString = "Feed", fontSize = 72), Text(extent = {{-38, 86}, {108, 2}}, lineColor = {255, 255, 255},
              lineThickness =                                                                                                                                                                                                        0.5, fillColor = {0, 0, 255},
              fillPattern =                                                                                                                                                                                                        FillPattern.Solid, fontSize = 72, textString = "Permeate"), Line(points = {{-100, 100}, {-100, -100}, {100, -100}, {100, 100}, {-100, 100}}, color = {0, 0, 0}, thickness = 0.5)}));
  end RO_middle_stream;

  model TestModel
    RO_hydraulics.Sources.Source_p_c source_p_c1(c = 400) annotation (
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Sources.Source_m_flow source_m_flow1(m_flow = -1 / 60) annotation (
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_hydraulics.TestConductor testConductor1 annotation (
      Placement(visible = true, transformation(origin = {2, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(testConductor1.port_b, source_m_flow1.port_a) annotation (
      Line(points = {{12, 0}, {90, 0}, {90, 0}, {90, 0}}, color = {0, 0, 255}));
    connect(source_p_c1.port_b, testConductor1.port_a) annotation (
      Line(points = {{-90, 0}, {-8, 0}, {-8, 0}, {-8, 0}}, color = {0, 0, 255}));
  end TestModel;

  model TestConductor
    extends RO_hydraulics.Interfaces.PartialTwoPort;
    parameter Real k(unit = "m.s") = 1e-5 "conductance in kg/(s.Pa)";
    Modelica.SIunits.MolarDensity port_a_c_inflow "concentration at port_a if inflowing";
    Modelica.SIunits.MolarDensity port_b_c_inflow "concentration at port_b if inflowing";
    parameter Modelica.SIunits.MolarDensity delta_c = 100 "concentration drop";
  equation
    // concentration for inflowing fluid
    port_a_c_inflow = inStream(port_a.c);
    port_b_c_inflow = inStream(port_b.c);
    // mass balance
    port_a.m_flow + port_b.m_flow = 0;
    // no change in concentration
    port_a.c = if m_flow > 0 then port_a_c_inflow else port_b_c_inflow - delta_c;
    port_b.c = if m_flow > 0 then port_a_c_inflow - delta_c else port_b_c_inflow;
    //  port_a.c*m_flow+port_b.c*m_flow = 0;
    // flow-pressure relationship
    m_flow = k * dp;
    annotation (
      Icon(graphics={  Rectangle(lineColor = {85, 0, 255}, fillColor = {85, 255, 255},
              fillPattern =                                                                          FillPattern.Solid, extent = {{-74, 22}, {74, -22}})}));
  end TestConductor;

  model RO_module "Series of segments of one RO module modeled as one element"
    Modelica.SIunits.EnergyDensity SpecificEnergy;
    RO_middle_stream rO_middle1(dx = 0.88 / 5) annotation (
      Placement(visible = true, transformation(origin = {-40, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle_stream rO_middle2(dx = 0.88 / 5) annotation (
      Placement(visible = true, transformation(origin = {0, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle_stream rO_middle3(dx = 0.88 / 5) annotation (
      Placement(visible = true, transformation(origin = {38, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle_stream rO_first(dx = 0.88 / 5) annotation (
      Placement(visible = true, transformation(origin = {-74, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    RO_middle_stream rO_last(dx = 0.88 / 5) annotation (
      Placement(visible = true, transformation(origin = {74, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Interfaces.Port_b permeate annotation (
      Placement(transformation(extent = {{90, -8}, {110, 12}}), iconTransformation(extent = {{90, -8}, {110, 12}})));
    Interfaces.Port_a feed annotation (
      Placement(transformation(extent = {{-112, -10}, {-92, 10}}), iconTransformation(extent = {{-112, -10}, {-92, 10}})));
    Interfaces.Port_b reject annotation (
      Placement(transformation(extent = {{-12, -110}, {8, -90}}), iconTransformation(extent = {{-12, -110}, {8, -90}})));
    Sources.Source_m_flow source_m_flow(m_flow = 0, c = 0) annotation (
      Placement(transformation(extent = {{-68, 56}, {-48, 76}})));
  equation
    SpecificEnergy = rO_first.port_feed_a.m_flow / rO_last.port_permeate_b.m_flow * rO_first.port_feed_b.p;
    connect(rO_first.port_feed_b, rO_middle1.port_feed_a) annotation (
      Line(points = {{-74, 0}, {-74, -10}, {-56, -10}, {-56, 10}, {-50, 10}}, color = {0, 0, 255}));
    connect(rO_first.port_permeate_b, rO_middle1.port_permeate_a) annotation (
      Line(points = {{-64, 10}, {-60, 10}, {-60, 24}, {-40, 24}, {-40, 20}}, color = {0, 0, 255}));
    connect(rO_middle1.port_feed_b, rO_middle2.port_feed_a) annotation (
      Line(points = {{-40, 0}, {-40, -10}, {-20, -10}, {-20, 10}, {-10, 10}}, color = {0, 0, 255}));
    connect(rO_middle2.port_feed_b, rO_middle3.port_feed_a) annotation (
      Line(points = {{0, 0}, {0, -8}, {20, -8}, {20, 10}, {28, 10}}, color = {0, 0, 255}));
    connect(rO_middle3.port_feed_b, rO_last.port_feed_a) annotation (
      Line(points = {{38, 0}, {38, -8}, {56, -8}, {56, 10}, {64, 10}}, color = {0, 0, 255}));
    connect(rO_middle1.port_permeate_b, rO_middle2.port_permeate_a) annotation (
      Line(points = {{-30, 10}, {-24, 10}, {-24, 26}, {0, 26}, {0, 20}}, color = {0, 0, 255}));
    connect(rO_middle2.port_permeate_b, rO_middle3.port_permeate_a) annotation (
      Line(points = {{10, 10}, {16, 10}, {16, 24}, {38, 24}, {38, 20}}, color = {0, 0, 255}));
    connect(rO_middle3.port_permeate_b, rO_last.port_permeate_a) annotation (
      Line(points = {{48, 10}, {52, 10}, {52, 28}, {74, 28}, {74, 20}}, color = {0, 0, 255}));
    connect(rO_first.port_feed_a, feed) annotation (
      Line(points = {{-84, 10}, {-102, 10}, {-102, 0}}, color = {0, 0, 255}));
    connect(rO_last.port_permeate_b, permeate) annotation (
      Line(points = {{84, 10}, {88, 10}, {88, 2}, {100, 2}}, color = {0, 0, 255}));
    connect(rO_last.port_feed_b, reject) annotation (
      Line(points = {{74, 0}, {74, -68}, {-2, -68}, {-2, -100}}, color = {0, 0, 255}));
    connect(source_m_flow.port_a, rO_first.port_permeate_a) annotation (
      Line(points = {{-68, 66}, {-74, 66}, {-74, 20}}, color = {0, 0, 255}));
    annotation (
      Icon(graphics={  Polygon(lineColor = {0, 0, 255}, fillColor = {7, 0, 255},
              fillPattern =                                                                    FillPattern.Solid, points = {{-100, 100}, {100, 100}, {100, -100}, {100, -100}, {-100, 100}}), Line(points = {{-100, -100}, {-100, 100}, {100, 100}, {100, -100}, {-100, -100}})}));
  end RO_module;

  model TestRO_module
    RO_hydraulics.Sources.Source_p_c source_permeate_p(c = 30) annotation (
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Sources.Source_m_flow source_permeate_flow(m_flow = 0 / 60) annotation (
      Placement(visible = true, transformation(origin = {-24, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    RO_middle_stream rO_middle_stream annotation (
      Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
    RO_middle_stream rO_middle_stream1 annotation (
      Placement(transformation(extent = {{-36, -10}, {-16, 10}})));
    RO_middle_stream rO_middle_stream2 annotation (
      Placement(transformation(extent = {{0, -10}, {20, 10}})));
    RO_middle_stream rO_middle_stream3 annotation (
      Placement(transformation(extent = {{30, -10}, {50, 10}})));
    RO_middle_stream rO_middle_stream4 annotation (
      Placement(transformation(extent = {{60, -10}, {80, 10}})));
  Modelica.Blocks.Sources.Cosine signal_pressure(                  freqHz = 2,
      amplitude=0,
      offset=35e5)                                                                             annotation (
      Placement(visible = true, transformation(origin={46,-78},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Cosine signal_mflow(                 freqHz = 2,
      amplitude=0.05,
      offset=0.08)                                                                         annotation (
      Placement(visible = true, transformation(origin={-134,0},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Sources.Source_p_c_signal source_reject_p annotation (Placement(visible=true,
          transformation(
          origin={80,-78},
          extent={{-10,-10},{10,10}},
          rotation=0)));
  RO_hydraulics.Sources.Source_m_flow_signal source_feed_m_flow annotation (
        Placement(visible=true, transformation(
          origin={-96,0},
          extent={{10,-10},{-10,10}},
          rotation=0)));
  equation
    connect(source_feed_m_flow.u, signal_mflow.y)
      annotation (Line(points={{-105.6,0},{-123,0}}, color={0,0,127}));
    connect(signal_pressure.y, source_reject_p.u)
      annotation (Line(points={{57,-78},{70.2,-78}}, color={0,0,127}));
    connect(rO_middle_stream.port_feed_b, rO_middle_stream1.port_feed_a) annotation (
      Line(points = {{-60, -10}, {-62, -10}, {-62, -24}, {-36, -24}, {-36, 0}}, color = {0, 0, 255}));
    connect(rO_middle_stream.port_permeate_b, rO_middle_stream1.port_permeate_a) annotation (
      Line(points = {{-50, 0}, {-42, 0}, {-42, 16}, {-26, 16}, {-26, 10}}, color = {0, 0, 255}));
    connect(rO_middle_stream.port_permeate_a, source_permeate_flow.port_a) annotation (
      Line(points = {{-60, 10}, {-60, 36}, {-14, 36}, {-14, 68}}, color = {0, 0, 255}));
    connect(rO_middle_stream1.port_feed_b, rO_middle_stream2.port_feed_a) annotation (
      Line(points = {{-26, -10}, {-24, -10}, {-24, -20}, {-2, -20}, {-2, -6}, {0, -6}, {0, 0}}, color = {0, 0, 255}));
    connect(rO_middle_stream1.port_permeate_b, rO_middle_stream2.port_permeate_a) annotation (
      Line(points = {{-16, 0}, {-6, 0}, {-6, 18}, {10, 18}, {10, 10}}, color = {0, 0, 255}));
    connect(rO_middle_stream2.port_feed_b, rO_middle_stream3.port_feed_a) annotation (
      Line(points = {{10, -10}, {10, -22}, {30, -22}, {30, 0}}, color = {0, 0, 255}));
    connect(rO_middle_stream2.port_permeate_b, rO_middle_stream3.port_permeate_a) annotation (
      Line(points = {{20, 0}, {24, 0}, {24, 16}, {40, 16}, {40, 10}}, color = {0, 0, 255}));
    connect(rO_middle_stream3.port_feed_b, rO_middle_stream4.port_feed_a) annotation (
      Line(points = {{40, -10}, {40, -22}, {56, -22}, {56, 0}, {60, 0}}, color = {0, 0, 255}));
    connect(rO_middle_stream3.port_permeate_b, rO_middle_stream4.port_permeate_a) annotation (
      Line(points = {{50, 0}, {52, 0}, {52, 20}, {70, 20}, {70, 10}}, color = {0, 0, 255}));
    connect(rO_middle_stream4.port_permeate_b, source_permeate_p.port_b) annotation (
      Line(points = {{80, 0}, {90, 0}}, color = {0, 0, 255}));
    connect(source_feed_m_flow.port_a, rO_middle_stream.port_feed_a)
      annotation (Line(points={{-86,0},{-70,0}}, color={0,0,255}));
    connect(source_reject_p.port_b, rO_middle_stream4.port_feed_b) annotation (
        Line(points={{90,-78},{96,-78},{96,-38},{70,-38},{70,-10}}, color={0,0,
            255}));
  end TestRO_module;

  package Interfaces
    connector Port_a "fluid pin for RO"
      Modelica.SIunits.Pressure p "pressure";
      stream Modelica.SIunits.MolarDensity c "Molar salt concentration";
      flow Modelica.SIunits.MassFlowRate m_flow "volumetric flow rate";
      annotation (
        defaultComponentName = "pin_n",
        Documentation(info = "<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>", revisions = "<html>
<dl>
<dt><i>1998</i></dt>
<dd>by Christoph Clauss initially implemented
</dd>
</dl>
</html>"),
        Icon(coordinateSystem(initialScale = 0.1), graphics={  Rectangle(lineColor = {0, 0, 255}, fillColor = {85, 0, 255},
                fillPattern =                                                                                                             FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                                                                                                                 FillPattern.Solid), Text(extent = {{-96, 100}, {104, 40}}, lineColor = {0, 0, 255}, textString = "%name")}));
    end Port_a;

    connector Port_b "fluid pin for RO"
      Modelica.SIunits.Pressure p "pressure";
      stream Modelica.SIunits.MolarDensity c "Molar salt concentration";
      flow Modelica.SIunits.MassFlowRate m_flow "volumetric flow rate";
      annotation (
        defaultComponentName = "pin_n",
        Documentation(info = "<html>
<p>Connectors PositivePin and NegativePin are nearly identical. The only difference is that the icons are different in order to identify more easily the pins of a component. Usually, connector PositivePin is used for the positive and connector NegativePin for the negative pin of an electrical component.</p>
</html>", revisions = "<html>
<dl>
<dt><i>1998</i></dt>
<dd>by Christoph Clauss initially implemented
</dd>
</dl>
</html>"),
        Icon(coordinateSystem(initialScale = 0.1), graphics={  Rectangle(lineColor = {0, 0, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                                FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}),
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                                                                                                                 FillPattern.Solid), Text(extent = {{-96, 100}, {104, 40}}, lineColor = {0, 0, 255}, textString = "%name")}));
    end Port_b;

    partial model PartialTwoPort
      Interfaces.Port_a port_a annotation (
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.Port_b port_b annotation (
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.SIunits.MassFlowRate m_flow "flow through element in positive direction";
      Modelica.SIunits.Pressure dp "pressure across element";
    protected
      parameter Boolean port_a_exposesState = false "= true if port_a exposes the state of a fluid volume";
      parameter Boolean port_b_exposesState = false "= true if port_b.p exposes the state of a fluid volume";
      parameter Boolean showDesignFlowDirection = true "= false to hide the arrow in the model icon";
    equation
      m_flow = port_a.m_flow;
      dp = port_a.p - port_b.p;
      annotation (
        Documentation(info = "<html>
  <p>
  This partial model defines an interface for components with two ports.
  The treatment of the design flow direction and of flow reversal are predefined based on the parameter <code><b>allowFlowReversal</b></code>.
  The component may transport fluid and may have internal storage for a given fluid <code><b>Medium</b></code>.
  </p>
  <p>
  An extending model providing direct access to internal storage of mass or energy through port_a or port_b
  should redefine the protected parameters <code><b>port_a_exposesState</b></code> and <code><b>port_b_exposesState</b></code> appropriately.
  This will be visualized at the port icons, in order to improve the understanding of fluid model diagrams.
  </p>
  </html>"),
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Polygon(points = {{20, -70}, {60, -85}, {20, -100}, {20, -70}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255},
                fillPattern =                                                                                                                                                                                                        FillPattern.Solid, visible = showDesignFlowDirection), Polygon(points = {{20, -75}, {50, -85}, {20, -95}, {20, -75}}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
                fillPattern =                                                                                                                                                                                                        FillPattern.Solid, visible = allowFlowReversal), Line(points = {{55, -85}, {-60, -85}}, color = {0, 128, 255}, visible = showDesignFlowDirection), Text(extent = {{-149, -114}, {151, -154}}, lineColor = {0, 0, 255}, textString = "%name"), Ellipse(extent = {{-110, 26}, {-90, -24}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                                                                                                                        FillPattern.Solid, visible = port_a_exposesState), Ellipse(extent = {{90, 25}, {110, -25}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0},
                fillPattern =                                                                                                                                                                                                        FillPattern.Solid, visible = port_b_exposesState)}));
    end PartialTwoPort;

    partial model PartialFourPort
      Interfaces.Port_a port_feed_a annotation (
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.Port_b port_feed_b annotation (
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.SIunits.MassFlowRate m_flow_reject "reject flow rate";
      Modelica.SIunits.MassFlowRate m_flow_permeate "permeate flow rate";
      Modelica.SIunits.MassFlowRate m_flow_membrane "flow rate across membrane";
      Modelica.SIunits.Pressure dp_feed "pressure difference across feed channel";
      Modelica.SIunits.Pressure dp_permeate "pressure difference across permeate channel";
    protected
      parameter Boolean showDesignFlowDirection = true "= false to hide the arrow in the model icon";
    public
      Interfaces.Port_a port_permeate_a annotation (
        Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.Port_b port_permeate_b annotation (
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      m_flow_reject = -port_feed_b.m_flow;
      m_flow_permeate = -port_permeate_b.m_flow;
      //mass balance on feed side
      m_flow_membrane = port_feed_a.m_flow + port_feed_b.m_flow;
      //mass balance on permeate side
      -m_flow_membrane = port_permeate_a.m_flow + port_permeate_b.m_flow;
      dp_feed = port_feed_a.p - port_feed_b.p;
      dp_permeate = port_permeate_a.p - port_permeate_b.p;
      annotation (
        Documentation(info = "<html>
  <p>
  This partial model defines an interface for components with four ports.
  The treatment of the design flow direction and of flow reversal are predefined based on the parameter <code><b>allowFlowReversal</b></code>.
  The component may transport fluid and may have internal storage for a given fluid <code><b>Medium</b></code>.
  </p>
  <p>
  An extending model providing direct access to internal storage of mass or energy through port_a or port_b
  should redefine the protected parameters <code><b>port_a_exposesState</b></code> and <code><b>port_b_exposesState</b></code> appropriately.
  This will be visualized at the port icons, in order to improve the understanding of fluid model diagrams.
  </p>
  </html>"),
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Polygon(points = {{-32, -54}, {-18, -84}, {-50, -76}, {-32, -54}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255},
                fillPattern =                                                                                                                                                                                                        FillPattern.Solid, visible = showDesignFlowDirection), Line(points = {{-32, -74}, {-86, -18}}, color = {0, 128, 255}, visible = showDesignFlowDirection), Text(extent = {{-149, -114}, {151, -154}}, lineColor = {0, 0, 255}, textString = "%name"), Line(points = {{74, 34}, {20, 90}}, color = {0, 128, 255}, visible = showDesignFlowDirection), Polygon(points = {{74, 54}, {88, 24}, {56, 32}, {74, 54}}, lineColor = {0, 128, 255}, fillColor = {0, 128, 255},
                fillPattern =                                                                                                                                                                                                        FillPattern.Solid, visible = showDesignFlowDirection)}));
    end PartialFourPort;
  end Interfaces;

  package Sources
    model Source_p_c
      Interfaces.Port_b port_b annotation (
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter Modelica.SIunits.Pressure p = 1E5 "pressure";
      parameter Modelica.SIunits.MolarDensity c = 600 "molar salt concentration, if flow goes out of element";
    equation
      port_b.p = p;
      port_b.c = if port_b.m_flow > 0 then inStream(port_b.c) else c;
      annotation (
        Icon(graphics={  Ellipse(origin = {0, -3}, fillColor = {85, 255, 255},
                fillPattern =                                                                FillPattern.Solid, extent = {{-58, 61}, {58, -57}}, endAngle = 360)}));
    end Source_p_c;

    model Source_m_flow
      parameter Modelica.SIunits.MassFlowRate m_flow = 1 / 60 "mass flow rate";
      parameter Modelica.SIunits.MolarDensity c = 600 "molar salt concentration";
      Interfaces.Port_a port_a annotation (
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      port_a.m_flow = -m_flow;
      port_a.c = if m_flow > 0 then c else inStream(port_a.c);
      annotation (
        Icon(graphics={  Ellipse(origin = {0, -3}, fillColor = {85, 255, 255},
                fillPattern =                                                                FillPattern.Solid, extent = {{-58, 61}, {58, -57}}, endAngle = 360), Line(origin = {-49.75, -82.25}, points = {{49.7481, 0.248069}, {-50.2519, 0.248069}, {-36.2519, 8.24807}, {-50.2519, 0.248069}, {-36.2519, -7.75193}}, color = {85, 0, 255})}, coordinateSystem(initialScale = 0.1)));
    end Source_m_flow;

    model Source_m_flow_signal
      parameter Modelica.SIunits.MolarDensity c = 600 "molar salt concentration";
      Interfaces.Port_a port_a annotation (
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput u annotation (
        Placement(transformation(extent = {{20, -20}, {-20, 20}}, rotation = 0, origin = {106, 0}), iconTransformation(extent = {{116, -20}, {76, 20}})));
    equation
      port_a.m_flow = -u;
      port_a.c = if u > 0 then c else inStream(port_a.c);
      annotation (
        Icon(graphics={  Ellipse(origin = {0, -3}, fillColor = {85, 255, 255},
                fillPattern =                                                                FillPattern.Solid, extent = {{-58, 61}, {58, -57}}, endAngle = 360), Line(origin = {-49.75, -82.25}, points = {{49.7481, 0.248069}, {-50.2519, 0.248069}, {-36.2519, 8.24807}, {-50.2519, 0.248069}, {-36.2519, -7.75193}}, color = {85, 0, 255})}, coordinateSystem(initialScale = 0.1)));
    end Source_m_flow_signal;


    model Source_p_c_signal
      Interfaces.Port_b port_b annotation (
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter Modelica.SIunits.MolarDensity c = 600 "molar salt concentration, if flow goes out of element";
      Modelica.Blocks.Interfaces.RealInput u annotation (
        Placement(visible = true, transformation(origin = {-98, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-98, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation
      port_b.p = u;
      port_b.c = if port_b.m_flow > 0 then inStream(port_b.c) else c;
      annotation (
        Icon(graphics={  Ellipse(origin = {0, -3}, fillColor = {85, 255, 255},
                fillPattern =                                                                FillPattern.Solid, extent = {{-58, 61}, {58, -57}}, endAngle = 360)}));
    end Source_p_c_signal;
  end Sources;
  annotation (
    uses(OpenHydraulics(version = "1.0"), Modelica(version = "3.2.2"), RO_middle(version = "1")));
end RO_hydraulics;
