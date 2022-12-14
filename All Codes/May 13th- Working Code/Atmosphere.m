function [T_atm, p_atm, rho, Mach] = Atmosphere(X)

global Velocity Altitude
if (Altitude<= 11000) 
        T_atm = 288.15-0.0065*Altitude;
        p_atm = 101325*(T_atm/288.15)^(9.81/(287*0.0065)); 
else 
        T_atm = 216.65;
        p_atm = 22632*exp(-9.81*(Altitude-11000)/(287*216.65)); 
end    
rho = p_atm/(287*T_atm);
    
Vsound=sqrt(1.4*p_atm/rho);
Mach=Velocity/Vsound;

end

