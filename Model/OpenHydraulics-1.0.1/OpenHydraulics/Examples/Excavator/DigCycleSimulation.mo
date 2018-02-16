within OpenHydraulics.Examples.Excavator;
model DigCycleSimulation
  extends Modelica.Icons.Example;
  OpenHydraulics.Examples.Excavator.SubSystems.MechanicsBody_noFriction mechanicsBody(swing_phi_start = 0, swing_initType = Modelica.Mechanics.MultiBody.Types.Init.PositionVelocity, boom_phi_start = 0.87266462599716, arm_phi_start = -1.3962634015955, bucket_phi_start = 0.34906585039887) annotation(Placement(visible = true, transformation(extent = {{10, -33}, {98, 46}}, rotation = 0)));
  OpenHydraulics.Examples.Excavator.SubSystems.HydraulicsSubSyst hydraulics(
      redeclare OpenHydraulics.Fluids.GenericOilSimple oil)
    annotation (Placement(transformation(extent={{-6,-18},{-46,22}},
          rotation=0)));
  inner Modelica.Mechanics.MultiBody.World world
    annotation (Placement(transformation(extent={{10,-60},{30,-40}},
          rotation=0)));
  OpenHydraulics.Examples.Excavator.SubSystems.DigCycleSeq digCycleSeq
    annotation (Placement(transformation(extent={{-98,-18},{-56,22}},
          rotation=0)));
equation
  connect(hydraulics.SwingFlange, mechanicsBody.swingFlange) annotation(Line(points = {{-9.2, 22}, {-9.2, 30}, {4, 30}, {4, -18}, {15, -18}}));
  connect(hydraulics.BucketCylRod, mechanicsBody.cylBucketRod) annotation(Line(points = {{-41.2, 22}, {-42, 22}, {-42, 58}, {86, 58}, {86, 34}}, color = {0, 127, 0}));
  connect(hydraulics.BucketCylBase, mechanicsBody.cylBucketBase) annotation(Line(points = {{-37.2, 22}, {-38, 22}, {-38, 54}, {82, 54}, {82, 34}}, color = {0, 127, 0}));
  connect(hydraulics.BoomCylBaseR, mechanicsBody.cylBoomRightBase) annotation(Line(points = {{-17.2, 22}, {-18, 22}, {-18, 34}, {0, 34}, {0, -30}, {46, -30}, {46, -21}, {51, -21}}, color = {0, 127, 0}));
  connect(hydraulics.BoomCylRodR, mechanicsBody.cylBoomRightRod) annotation(Line(points = {{-21.2, 22}, {-21.2, 38}, {2, 38}, {2, -26}, {55, -26}, {55, -21}}, color = {0, 127, 0}));
  connect(hydraulics.BoomCylBaseL, mechanicsBody.cylBoomLeftBase) annotation(Line(points = {{-6, 10.8}, {15, 10.8}, {15, 12}}, color = {0, 127, 0}));
  connect(hydraulics.BoomCylRodL, mechanicsBody.cylBoomLeftRod) annotation(Line(points = {{-6, 14.8}, {15, 14.8}, {15, 16}}, color = {0, 127, 0}));
  connect(world.frame_b, mechanicsBody.baseFrame) annotation(Line(points = {{30, -50}, {40, -50}, {40, -20}}, color = {95, 95, 95}, thickness = 0.5));
  connect(hydraulics.ArmCylBase, mechanicsBody.cylArmBase) annotation(Line(points = {{-26.8, 22}, {-26, 22}, {-26, 46}, {44, 46}, {44, 34}, {49, 34}}, color = {0, 127, 0}));
  connect(hydraulics.ArmCylRod, mechanicsBody.cylArmRod) annotation(Line(points = {{-30.8, 22}, {-30, 22}, {-30, 50}, {48, 50}, {48, 34}, {53, 34}}, color = {0, 127, 0}));
  connect(digCycleSeq.y1, hydraulics.Command) annotation(Line(points = {{-53.9, 2}, {-50.85, 2}, {-47.8, 2}}, color = {0, 0, 127}));
  annotation (Diagram(graphics),
    experiment(StopTime=20, Tolerance=1e-008),
    experimentSetupOutput);
end DigCycleSimulation;