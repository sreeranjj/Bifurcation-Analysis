%===============%
% Trim Equation %
%===============%
function F = Trim(x)
d2r = pi/180;   % deg to rad conversion
r2d = 180/pi;   % rad to deg conversion

global Velocity Altitude 

% Trim - Independent variables are
% x(1) = alpha
% x(2) = beta
% x(3) = del_e
% x(4) = del_a
% x(5) = del_r
% x(6) = Thrtl

symetri = 0;
% For  Wing level trim Lateral Directional Independent variables are made 0
if (symetric == 1)
    x(2) = 0.0;
    x(4) = 0.0;
    x(5) = 0.0;
end


% Reading Aero and Geometric data
[Geom, ALPHA_BREAK, F18_Aerodata] = Aerodata;

State   =  [Velocity; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; Altitude];  % Initial Conditions

[T_atm, p_atm, rho, Mach, g] = Atmosphere(State);    % Atmosphere model

%% Function
%=========================================================================%
%         Independent variables Assignment to local variables             %
%=========================================================================%
Alpha = x(1);
Beta  = x(2);
Del_e = x(3);
Del_a = x(4);
Del_r = x(5);
Thrtl = x(6);

% In radians
Alphar = Alpha*d2r;
Betar  = Beta*d2r;
Del_er = Del_e*d2r;
Del_ar = Del_a*d2r; 
Del_rr = Del_r*d2r;

Theta = Alpha;
Phi   = 0.0;
Psi   = 0.0;

State =  [Velocity; Alphar; Betar; p; q; r; Phir; Thetar; Psir; 0; 0; Altitude];  
cntl  =  [Del_er; Del_er; Del_ar; Del_rr; Thrtl; Thrtl; 0; 0; 0; 0];

[T_atm, p_atm, rho, Mach, g] = Atmosphere(State);   % Atmosphere model

[Thrust] = Engine(0.0, 1, State, cntl, Mach, g);    % Engine model

[F] = Equations_of_Motion(State,g,ALPHA_BREAK,F18_Aerodata,Thrust,Geom,rho,cntl); % Aircraft model

Fun = [F(1); F(2)*r2d; F(5)*r2d; F(3)*r2d; F(4)*r2d; F(6)*r2d; ...
       F(7)*r2d; F(8)*r2d; F(9)*r2d;  F(10); F(11); F(12)];

Ftrim = Fun(1:6);     % Xdot in trim equation

F = Ftrim;
