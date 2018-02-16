model RO_first
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
  rosmosis=A_w*(p_fin - p_pout) ;
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
if (p_fin - delta_p)>0
 then p_fout = p_fin - delta_p;
  else
   p_fout=0;
   end if;
//Power=Q_fin*p_fin/Q_pout /(36); //kWh/m3 of desalted water
  annotation(uses(Modelica(version="3.2.2")), Diagram(graphics={
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