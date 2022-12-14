function [Geom] = Geom_data()

m = 36100/2.2046; % 36100 lbs
mo=50;
m_tot=m+mo;
g = 9.81; % m/s5
W = m*g; % N
Ixx = 22789*14.59/10.764; % slugs ft2 -> Kg m2
Iyy = 176809*14.59/10.764;
Izz = 191744*14.59/10.764;
Ixz = -2305*14.59/10.764;
Ixy=-305*14.59/10.764;
Iyz=-505*14.59/10.764;
rho_sl = 1.225; % density @ sea-level
S = 400.0*0.3048*0.3048; % ft2 -> m2 (Wing planform area)
b = 37.4*0.3048; % ft -> m (Wing Span)
cbar = 11.52*0.3048; % ft -> m; Mean Aerodynamic Chord (MAC)
max_thrust = 32000*g/2.2046; % F-18 max. thrust (sea-level) static: 16000 lbs x 2
xcg = 0.238*cbar; % CG location referenced from the leading edge of the mean aerodynamic chord 
xcm=0.05;
ycm=0.1;
zcm=0;

Geom=zeros(20,1);

Geom(1,1)=m;
Geom(2,1)=mo;
Geom(3,1)=m_tot;
Geom(4,1)=g;
Geom(5,1)=W;
Geom(6,1)=Ixx;
Geom(7,1)=Iyy;
Geom(8,1)=Izz;
Geom(9,1)=Ixy;
Geom(10,1)=Ixz;
Geom(11,1)=Iyz;
Geom(12,1)=rho_sl;
Geom(13,1)=S;
Geom(14,1)=b;
Geom(15,1)=cbar;
Geom(16,1)=max_thrust;
Geom(17,1)=xcg;
Geom(18,1)=xcm;
Geom(19,1)=ycm;
Geom(20,1)=zcm;

end




