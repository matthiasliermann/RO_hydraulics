model RO

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
  connect(rO_middle3.pin_n5, rO_last1.pin_n1) annotation(Line(points = {{60, 10}, {68, 10}, {68, 28}, {90, 28}, {90, 20}, {90, 20}, {90, 20}}, color = {0, 0, 255}));
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
  
  annotation(uses(RO_middle(version = "1"), OpenHydraulics(version = "1.0")), Icon(graphics = {Polygon(lineColor = {0, 0, 255}, fillColor = {7, 0, 255}, fillPattern = FillPattern.Solid, points = {{-100, 100}, {100, 100}, {100, -100}, {100, -100}, {-100, 100}}), Line(points = {{-100, -100}, {-100, 100}, {100, 100}, {100, -100}, {-100, -100}})}), Diagram);
end RO;