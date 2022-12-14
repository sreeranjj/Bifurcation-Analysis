clear all
% clc

global ALPHA_BREAK clift clift_de cd0 cd_de cm cm_de cy_de croll_de cn_de m g S rho Tm V gamma Ixx Iyy Izz Ixz b cbar

[ALPHA_BREAK,clift,clift_q,cd0,cd_q,cd_de,cy_b,croll_b,croll_p,croll_r,cm,cm_q,cn_b,cn_p,cy_r,...
    clift_de,cy_da,cy_dr,croll_da,cn_da,cn_r,cy_p,croll_dr,cn_dr,cm_de,cy_de,croll_de,cn_de] = Aerodata1;


% Geometric / mass properties data
m = 36100/2.2046; % 36100 lbs
g = 9.81; % m/s2
W = m*g; % N
rho_sl = 1.225; % density @ sea-level
S = 400.0*0.3048*0.3048; % ft2 -> m2 (Wing planform area)
max_thrust = 32000*g/2.2046; % F-18 max. thrust (sea-level) static: 16000 lbs x 2
Ixx = 22789*14.59/10.764; % slugs ft2 -> Kg m2
Iyy = 176809*14.59/10.764;
Izz = 191744*14.59/10.764;
Ixz = -2305*14.59/10.764;
b = 37.4*0.3048; % ft -> m (Wing Span)
cbar = 11.52*0.3048; % ft -> m; Mean Aerodynamic Chord (MAC)

ALT=5000;%4500;   %enter altitude in m
V=120;%59.2+1.6;      %enter velocity in m/s
gamma=0;    % enter flight path angle in deg


if (ALT <= 11000) 
   T_atm = 288.15-0.0065*ALT;
   p_atm = 101325*(T_atm/288.15)^(9.81/(287*0.0065)); 
else 
   T_atm = 216.65;
   p_atm = 22632*exp(-9.81*(ALT-11000)/(287*216.65)); 
end    
rho = p_atm/(287*T_atm);

% M=0.6;
% Tm=0.9*max_thrust*(0.88+0.24*(M-0.6)^1.4)*(rho/rho_sl)^0.8;
Tm=((rho/rho_sl)^0.7)*max_thrust*(1-exp((ALT-17000)/2000)); %from Stengel
% Tm=max_thrust at the given altitude;

[x_trim,fval]=fsolve(@(x) Trim(x),[0.026;0;0;0;0;0.1])


