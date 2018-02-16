within OpenHydraulics.Components.MotorsPumps;
model VariableDisplacementPump "Variable Displacement Pump with losses"

  // parameters in sizing tab
  parameter SI.Volume Dmax = 0.001 "Maximum pump displacement"
    annotation(Dialog(tab="Sizing"));
  parameter SI.Volume Dmin = 0 "Minimum pump displacement (<0 for over-center)"
    annotation(Dialog(tab="Sizing"));
  parameter SI.Volume Dlimit = max(abs(Dmax),abs(Dmin)) "Displacement of pump"
    annotation(Dialog(tab="Sizing"));

  // loss parameters
  parameter Real Cv = 60000 "Coefficient of viscous drag"
    annotation(Dialog(tab="Efficiency",group="Mechanical Efficiency"));
  parameter Real Cf = 0.007
    "Coefficient of Coulomb friction (fraction of full stroke torque)"
    annotation(Dialog(tab="Efficiency",group="Mechanical Efficiency"));
  parameter Real Cs = 1.8e-9 "Leakage coefficient"
    annotation(Dialog(tab="Efficiency",group="Volumetric Efficiency"));
  parameter Real Vr = 0.54 "Volume ratio of pump or motor"
    annotation(Dialog(tab="Efficiency",group="Volumetric Efficiency"));
  parameter SI.DynamicViscosity mu = oil.dynamicViscosity(p_init)
    "Dynamic viscosity (used only for efficiency computation)"
    annotation(Dialog(tab="Efficiency",group="Related Parameters"));
  parameter SI.BulkModulus B = oil.approxBulkModulus(p_init)
    "Approximate bulk modulus (used only for efficiency models)"
    annotation(Dialog(tab="Efficiency",group="Related Parameters"));

  // friction
  BaseClasses.MechanicalPumpLosses mechanicalPumpLosses(
    Cv=Cv,
    Cf=Cf,
    Dmax=Dlimit,
    mu = mu)
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}},
          rotation=0)));

  // the components
  Basic.FluidPower2MechRotVar fluidPower2MechRot(
    final Dmax=Dmax,
    final Dmin=Dmin,
    final Dlimit=Dlimit)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=0)));
  Basic.VariableFlow leakage_P2T
    annotation (Placement(transformation(
        origin={40,0},
        extent={{-10,10},{10,-10}},
        rotation=270)));
  Lines.NJunction j1
    annotation (Placement(transformation(extent={{-10,-50},{10,-30}},
          rotation=0)));
  Lines.NJunction j2
    annotation (Placement(transformation(extent={{-10,30},{10,50}},
          rotation=0)));

  // the connectors
  OpenHydraulics.Interfaces.FluidPort portP
    annotation (Placement(transformation(extent={{-10,90},{10,110}},
          rotation=0)));
  OpenHydraulics.Interfaces.FluidPort portT
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}},
          rotation=0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a
    "(left) driving flange (flange axis directed INTO cut plane)"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealInput dispFraction
    annotation (Placement(transformation(extent={{-100,-96},{-68,-64}},
          rotation=0)));

  extends OpenHydraulics.Interfaces.PartialFluidComponent;
protected
  SI.Pressure dp = portP.p - portT.p;
equation
  // connect the input of the leakage model and the mechanical loss model
  leakage_P2T.q_flow = Cs*dp*Dlimit/mu
      + dp*mechanicalPumpLosses.w*Dlimit*(Vr + (1+dispFraction)/2)/B;
  mechanicalPumpLosses.dp = dp;

  connect(fluidPower2MechRot.dispFraction, dispFraction) annotation (Line(
        points={{-8.5,-7.9},{-42,-7.9},{-42,-80},{-84,-80}}, color={0,0,127}));
  connect(portT, j1.port[1]) annotation (Line(points={{0,-100},{0,-40.6667}},
        color={255,0,0}));
  connect(portP, j2.port[1])
    annotation (Line(points={{0,100},{0,39.3333}}, color={255,0,0}));
  connect(fluidPower2MechRot.port_b, j2.port[2])
    annotation (Line(points={{0,10},{0,40}}, color={255,0,0}));
  connect(fluidPower2MechRot.port_a, j1.port[2]) annotation (Line(points={{
          0,-10},{0,-40}}, color={255,0,0}));
  connect(leakage_P2T.port_a, j2.port[3]) annotation (Line(points={{40,10},
          {40,40},{0,40},{0,40.6667}}, color={255,0,0}));
  connect(leakage_P2T.port_b, j1.port[3]) annotation (Line(points={{40,-10},
          {40,-40},{0,-40},{0,-39.3333}}, color={255,0,0}));
  connect(flange_a, mechanicalPumpLosses.flange_a)
    annotation (Line(points={{-100,0},{-80,0}}, color={0,0,0}));
  connect(mechanicalPumpLosses.flange_b, fluidPower2MechRot.flange_a)
    annotation (Line(points={{-60,0},{-10,0}}, color={0,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
        Ellipse(
          extent={{-54,54},{54,-54}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{0,-54},{0,-100}}, color={255,0,0}),
        Line(points={{0,100},{0,54}}, color={255,0,0}),
        Rectangle(
          extent={{-90,8},{-54,-8}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(extent={{100,-54},{-100,-90}}, textString="%name"),
        Text(
          extent={{10,-80},{40,-120}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="T"),
        Polygon(
          points={{-20,34},{0,54},{20,34},{-20,34}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-20,-34},{0,-54},{20,-34},{-20,-34}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{10,120},{40,80}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="P"),
        Polygon(
          points={{80,80},{52,66},{66,52},{80,80}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{80,80}}, color={0,0,0})}),
                            Diagram(graphics={
        Text(
          extent={{52,76},{52,64}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString=
               "section for"),
        Text(
          extent={{52,88},{52,76}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString=
               "see equations"),
        Rectangle(extent={{12,88},{90,52}}, lineColor={0,0,255}),
        Text(
          extent={{52,64},{52,52}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString=
               "loss relationships")}));
end VariableDisplacementPump;
