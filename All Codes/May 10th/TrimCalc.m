

global Velocity Altitude
Velocity=120;
Altitude=5000;
[x_trim,fval]=fsolve(@(x) Trim(x),[0.026;0;0;0;0;0.1]);
x_trim
%options=optimoptions(@fsolve,'MaxFunEvals',10000,'Maxiter',2000);