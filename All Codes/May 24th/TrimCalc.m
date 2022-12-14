

global Velocity Altitude
Velocity=120;
Altitude=5000;

[x_trim,fval]=fsolve(@(x) Trim(x),[15;0.026;0;0;0;0.1],options)

options=optimoptions(@fsolve,'TolFun',1E-12,'TolX',1E-12);