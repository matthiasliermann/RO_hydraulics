model RO_system

//extends OpenHydraulics.Interfaces.PartialFluidComponent;

  extends OpenHydraulics.Interfaces.PartialFluidCircuit(redeclare
     OpenHydraulics.Fluids.GenericOilSimple oil);
  RO rO1 annotation(Placement(visible = true, transformation(origin = {-22, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  OpenHydraulics.Basic.VarPressureSource source annotation(Placement(visible = true, transformation(origin = {-58, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Cosine cosine1(amplitude = 35e5, freqHz = 2, offset = 36e5)  annotation(Placement(visible = true, transformation(origin = {-118, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  OpenHydraulics.Basic.VarVolumeSource source2 annotation(Placement(visible = true, transformation(origin = {-22, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Cosine cosine2(amplitude = 0e-5, freqHz = 2, offset = -20e-5) annotation(Placement(visible = true, transformation(origin = {-68, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(cosine2.y, source2.control) annotation(Line(points = {{-56, -28}, {-32, -28}, {-32, -28}, {-32, -28}}, color = {0, 0, 127}));
  connect(rO1.port1, source2.port) annotation(Line(points = {{-22, 8}, {-22, 8}, {-22, -18}, {-22, -18}}, color = {255, 0, 0}));
  connect(cosine1.y, source.control) annotation(Line(points = {{-106, 8}, {-68, 8}, {-68, 8}, {-68, 8}}, color = {0, 0, 127}));
  connect(source.port, rO1.port) annotation(Line(points = {{-58, 18}, {-32, 18}}, color = {255, 0, 0}));
  annotation(uses(RO_middle(version = "1"), OpenHydraulics(version = "1.0"), Modelica(version = "3.2.1")));
end RO_system;