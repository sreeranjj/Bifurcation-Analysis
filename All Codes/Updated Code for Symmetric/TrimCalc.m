

global Velocity Altitude
Velocity=100;
Altitude=1500;
options=optimoptions(@fsolve,'TolFun',1E-14,'TolX',1E-14);
[x_trim,fval]=fsolve(@(x) Trim(x),[15;0.026;0;0;0;0.1],options)

