function [Geom] = Geom()



m = 36100/2.2046; % 36100 lbs
mo=0;
m_tot=m+mo;
g = 9.81; % m/s5
W = m*g; % N
Ixx = 22789*14.59/10.764; % slugs ft2 -> Kg m2
Iyy = 176809*14.59/10.764;
Izz = 191744*14.59/10.764;
Ixz = -2305*14.59/10.764;
Ixy=0;
Iyz=0;
rho_sl = 1.225; % density @ sea-level
S = 400.0*0.3048*0.3048; % ft2 -> m2 (Wing planform area)
b = 37.4*0.3048; % ft -> m (Wing Span)
cbar = 11.52*0.3048; % ft -> m; Mean Aerodynamic Chord (MAC)
max_thrust = 32000*g/2.2046; % F-18 max. thrust (sea-level) static: 16000 lbs x 2
xcg = 0.238*cbar; % CG location referenced from the leading edge of the mean aerodynamic chord 

Geom=zeros(17);

Geom(1)=m;
Geom(2)=mo;
Geom(3)=m_tot;
Geom(4)=g;
Geom(5)=W;
Geom(6)=Ixx;
Geom(7)=Iyy;
Geom(8)=Izz;
Geom(9)=Ixy;
Geom(10)=Ixz;
Geom(11)=Iyz;
Geom(12)=rho_sl;
Geom(13)=S;
Geom(14)=b;
Geom(15)=cbar;
Geom(16)=max_thrust;
Geom(17)=xcg;




