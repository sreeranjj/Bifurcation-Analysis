% Propogation Starts here



d2r=pi/180; %deg to rad conversion
r2d=180/pi; %rad to deg conversion
global Velocity Altitude


symmetric=2;
Velocity=120;
if symmetric==1
alpha=8.497995890268339*d2r;%8.626914799287016*d2r;%;%5.0142;  %deg value converted to rad
beta=-0.000000054837134*d2r;%0.497561938330999*d2r; %deg value converted to rad
p=0*d2r;    %deg/sec value converted to rad/sec
q=0*d2r;    %deg/sec value converted to rad/sec
r=0*d2r;    %deg/sec value converted to rad/sec
phi=0*d2r;  %deg value converted to rad
theta=8.497995890268339*d2r;%8.626914799287016*d2r;%7.0142;   % theta=alpha+gamma  %deg value converted to rad
psi=0*d2r;  %deg value converted to rad
xcg = 0.238*cbar; % CG location referenced from the leading edge of the mean aerodynamic chord 
% mlgear_loc = 0.39*cbar; % Main landing gear location 
del_e=-0.916369995782296*d2r;%-1.655500079239364*d2r;%-0.2409;     %deg value converted to rad
del_a=-0.000000040270334*d2r;%-5.618532353068066*d2r;        %deg value converted to rad
del_r=-0.000000281141234*d2r;%1.888977775840741*d2r;        %deg value converted to rad
eta= 0.228927295626050;%0.234991293002048;%0.3404;
xe=0;
ye=0;
Altitude=5000;
end

if symmetric==0
Velocity=120;
alpha=8.626914799287016*d2r;%;%5.0142;  %deg value converted to rad
beta=0.497561938330999*d2r; %deg value converted to rad
p=0*d2r;    %deg/sec value converted to rad/sec
q=0*d2r;    %deg/sec value converted to rad/sec
r=0*d2r;    %deg/sec value converted to rad/sec
phi=0*d2r;  %deg value converted to rad
theta=8.626914799287016*d2r;%7.0142;   % theta=alpha+gamma  %deg value converted to rad
psi=0*d2r;  %deg value converted to rad
xcg = 0.238*cbar; % CG location referenced from the leading edge of the mean aerodynamic chord 
% mlgear_loc = 0.39*cbar; % Main landing gear location 
del_e=-1.655500079239364*d2r;%-0.2409;     %deg value converted to rad
del_a=-5.618532353068066*d2r;        %deg value converted to rad
del_r=1.888977775840741*d2r;        %deg value converted to rad
eta= 0.234991293002048;%0.3404;
xe=0;
ye=0;
Altitude=5000;
end

if symmetric==2
Velocity=1.0023427285E+02;
alpha=2.3274039490E-01;%;%5.0142;  %deg value converted to rad
beta=3.5050009547E-02; %deg value converted to rad
p=-4.6716399840E-03;    %deg/sec value converted to rad/sec
q=8.2274528671E-03;    %deg/sec value converted to rad/sec
r=2.6249527892E-02;    %deg/sec value converted to rad/sec
phi=3.0373410084E-01;  %deg value converted to rad
theta=1.6821917154E-01;%7.0142;   % theta=alpha+gamma  %deg value converted to rad
psi=0*d2r;  %deg value converted to rad
xcg = 0.238*cbar; % CG location referenced from the leading edge of the mean aerodynamic chord 
% mlgear_loc = 0.39*cbar; % Main landing gear location 
del_e=-3.3993921228E-02;%-0.2409;     %deg value converted to rad
del_a=-5.618532353068066*d2r;        %deg value converted to rad
del_r=1.888977775840741*d2r;        %deg value converted to rad
eta= 0.234991293002048;%0.3404;
xe=0;
ye=0;
Altitude=5000;
end

u= Velocity*cos(alpha)*cos(beta);
v= Velocity*sin(beta);
w= Velocity*sin(alpha)*cos(beta);

%IC=[V;alpha;beta;p;q;r;phi;theta;psi;xe;ye;ALT]; 
IC=[u;v;w;p;q;r;phi;theta;psi;xe;ye;Altitude];


cntl=[del_e; del_e; del_a; del_r; eta; eta; 0; 0; 0; 0];
t0=0;
tf=500;
Ts=0.1;  %step size=0.1sec
Geom=Geom_data;
x=IC;
X=zeros(1+(tf-t0)/Ts,length(IC));
X(1,:)=IC';
    Vtemp=sqrt(X(1,1)^2+X(1,2)^2+X(1,3)^2);
    alphatemp=atan(X(1,3)/X(1,1));
    betatemp=asin(X(1,2)/Vtemp);
    X(1,1)=Vtemp;
    X(1,2)=alphatemp;
    X(1,3)=betatemp;
ETA=zeros(tf/Ts,1);
tic;
for i=1:tf/Ts
   
    
    
    
%     if i>=100 
%         del_e=(-5.2406+2.0)*d2r;
%     else
%         del_e=0*d2r;
%     end

%     if i>=100 & i<120
%         eta=(0.6156+0.1);        
%     else
%         if i>=120 & i<140
%             eta=(0.6156-0.1);
%         else
%             eta=0.6156;
%         end         
%     end
    
%         if i>=100 
%             eta=(0.6156+0.1);
%         end



%     if i>=100 && i<120
%         del_r=3.0*d2r;
%     else
%         del_r=0;
%     end

%     if i>=100 && i<120
%         del_a=5*d2r;        
% %     else
% %         if i>=1000 && i<1020
% %             del_a=-5*d2r;  
% %         else
% %             del_a=0;
% %         end         
%     end
    
    ETA(i,1)=eta;
    
        
    if (x(12)<= 11000) 
        T_atm = 288.15-0.0065*x(12);
        p_atm = 101325*(T_atm/288.15)^(9.81/(287*0.0065)); 
    else 
        T_atm = 216.65;
        p_atm = 22632*exp(-9.81*(x(12)-11000)/(287*216.65)); 
    end    
    rho = p_atm/(287*T_atm);

    % Tm=0.9*max_thrust*(0.88+0.24*(M-0.6)^1.4)*(rho/rho_sl)^0.8;
    Thrust=((rho/rho_sl)^0.7)*max_thrust*(1-exp((x(12)-17000)/2000)); %from Stengel
    % Tm=max_thrust;
    

%     CD0=interp1(ALPHA_BREAK,cd0,x(2));
%     CD_q=interp1(ALPHA_BREAK,cd_q,x(2));
%     CD_de=interp1(ALPHA_BREAK,cd_de,x(2));
% 
%     CL0=interp1(ALPHA_BREAK,clift,x(2));
%     CL_q=interp1(ALPHA_BREAK,clift_q,x(2));
%     CL_de=interp1(ALPHA_BREAK,clift_de,x(2));
% 
%     CY_b=interp1(ALPHA_BREAK,cy_b,x(2)); 
%     CY_p=interp1(ALPHA_BREAK,cy_p,x(2)); 
%     CY_r=interp1(ALPHA_BREAK,cy_r,x(2)); 
%     CY_de=interp1(ALPHA_BREAK,cy_de,x(2)); 
%     CY_da=interp1(ALPHA_BREAK,cy_da,x(2)); 
%     CY_dr=interp1(ALPHA_BREAK,cy_dr,x(2));
% 
%     Cl_b=interp1(ALPHA_BREAK,croll_b,x(2)); 
%     Cl_p=interp1(ALPHA_BREAK,croll_p,x(2)); 
%     Cl_r=interp1(ALPHA_BREAK,croll_r,x(2)); 
%     Cl_de=interp1(ALPHA_BREAK,croll_de,x(2)); 
%     Cl_da=interp1(ALPHA_BREAK,croll_da,x(2)); 
%     Cl_dr=interp1(ALPHA_BREAK,croll_dr,x(2));
% 
%     Cn_b=interp1(ALPHA_BREAK,cn_b,x(2)); 
%     Cn_p=interp1(ALPHA_BREAK,cn_p,x(2)); 
%     Cn_r=interp1(ALPHA_BREAK,cn_r,x(2)); 
%     Cn_de=interp1(ALPHA_BREAK,cn_de,x(2)); 
%     Cn_da=interp1(ALPHA_BREAK,cn_da,x(2)); 
%     Cn_dr=interp1(ALPHA_BREAK,cn_dr,x(2));
% 
%     Cm0=interp1(ALPHA_BREAK,cm,x(2));
%     Cm_q=interp1(ALPHA_BREAK,cm_q,x(2));
%     Cm_de=interp1(ALPHA_BREAK,cm_de,x(2));
    
    
%     CD0=interp1(ALPHA_BREAK,cd0,x(2),'spline');
%     CD_q=interp1(ALPHA_BREAK,cd_q,x(2),'spline');
%     CD_de=interp1(ALPHA_BREAK,cd_de,x(2),'spline');
% 
%     CL0=interp1(ALPHA_BREAK,clift,x(2),'spline');
%     CL_q=interp1(ALPHA_BREAK,clift_q,x(2),'spline');
%     CL_de=interp1(ALPHA_BREAK,clift_de,x(2),'spline');
% 
%     CY_b=interp1(ALPHA_BREAK,cy_b,x(2),'spline'); 
%     CY_p=interp1(ALPHA_BREAK,cy_p,x(2),'spline'); 
%     CY_r=interp1(ALPHA_BREAK,cy_r,x(2),'spline'); 
%     CY_de=interp1(ALPHA_BREAK,cy_de,x(2),'spline'); 
%     CY_da=interp1(ALPHA_BREAK,cy_da,x(2),'spline'); 
%     CY_dr=interp1(ALPHA_BREAK,cy_dr,x(2),'spline');
% 
%     Cl_b=interp1(ALPHA_BREAK,croll_b,x(2),'spline'); 
%     Cl_p=interp1(ALPHA_BREAK,croll_p,x(2),'spline'); 
%     Cl_r=interp1(ALPHA_BREAK,croll_r,x(2),'spline'); 
%     Cl_de=interp1(ALPHA_BREAK,croll_de,x(2),'spline'); 
%     Cl_da=interp1(ALPHA_BREAK,croll_da,x(2),'spline'); 
%     Cl_dr=interp1(ALPHA_BREAK,croll_dr,x(2),'spline');
% 
%     Cn_b=interp1(ALPHA_BREAK,cn_b,x(2),'spline'); 
%     Cn_p=interp1(ALPHA_BREAK,cn_p,x(2),'spline'); 
%     Cn_r=interp1(ALPHA_BREAK,cn_r,x(2),'spline'); 
%     Cn_de=interp1(ALPHA_BREAK,cn_de,x(2),'spline'); 
%     Cn_da=interp1(ALPHA_BREAK,cn_da,x(2),'spline'); 
%     Cn_dr=interp1(ALPHA_BREAK,cn_dr,x(2),'spline');
% 
%     Cm0=interp1(ALPHA_BREAK,cm,x(2),'spline');
%     Cm_q=interp1(ALPHA_BREAK,cm_q,x(2),'spline');
%     Cm_de=interp1(ALPHA_BREAK,cm_de,x(2),'spline');

    x;
    k1=Ts.*Equations_of_Motion(t0+(i-1)*Ts,x,g,Thrust,Geom,rho,cntl);
    k2=Ts.*Equations_of_Motion((t0+(i-1)*Ts+Ts/2),(x+k1/2),g,Thrust,Geom,rho,cntl);
    k3=Ts.*Equations_of_Motion((t0+(i-1)*Ts+Ts/2),(x+k2/2),g,Thrust,Geom,rho,cntl);
    k4=Ts.*Equations_of_Motion((t0+(i-1)*Ts+Ts),(x+k3),g,Thrust,Geom,rho,cntl);
    x=x+(k1+2.*k2+2.*k3+k4)/6;
    
%     x=RK_4_rad((t0+(i-1)*Ts),Ts,x);
    
  
   
%     X(i+1,1)=x(1,1);
    for j=1:12
        X(i+1,j)=x(j,1) ;
    end
   
    
    Vtemp2=sqrt(X(i+1,1)^2+X(i+1,2)^2+X(i+1,3)^2);
    alphatemp2=atan(X(i+1,3)/X(i+1,1));
    betatemp2=asin(X(i+1,2)/Vtemp2);
    X(i+1,1)=Vtemp2;
    X(i+1,2)=alphatemp2;
    X(i+1,3)=betatemp2;
    
end
toc;

for j=2:9
    X(:,j)=X(:,j).*r2d;     %all angular variables are converted back to deg
end

t=t0:Ts:tf;



 
figure(1)
 for j=1:3
     subplot(3,1,j)
     plot(t,X(:,j))
     if j==1
         title('V (m/s) v/s Time(s)');
     end
     if j==2
         title('Alpha (degree) v/s Time(s)');
     end
     if j==3
         title('Beta (degree) v/s Time(s)');
     end
%      if IC(j,1)~=0
%         axis([t0 tf 0.95*IC(j,1) 1.05*IC(j,1)]);
%      else
%         axis([t0 tf -0.1 0.1]); 
%      end
     grid on
 end

figure(2)
 for j=1:3
     subplot(3,1,j)
     plot(t,X(:,j+3))
            if j==1
            title('p (deg/s) v/s Time(s)');
            end
            if j==2
            title('q (deg/s) v/s Time(s)');
            end
            if j==3
            title('r (deg/s) v/s Time(s)');
            end
%      if IC(j+3,1)~=0
%         axis([t0 tf 0.95*IC(j+3,1) 1.05*IC(j+3,1)]);
%      else
%         axis([t0 tf -0.1 0.1]);
%      end
     grid on
 end
 
figure(3)
 for j=1:3
     subplot(3,1,j)
     plot(t,X(:,j+6))
            if j==1
            title('phi (deg) v/s Time(s)');
            end
            if j==2
            title('theta (deg) v/s Time(s)');
            if symmetric==1
            ylim([8.49,8.50]);
            end
            end
            if j==3
            title('psi (deg) v/s Time(s)');
            end
%      if IC(j+6,1)~=0
%         axis([t0 tf 0.95*IC(j+6,1) 1.05*IC(j+6,1)]);
%      else
%         axis([t0 tf -0.1 0.1]);
%      end
     grid on
 end 

figure(4)
 for j=1:3
     subplot(3,1,j)
     plot(t,X(:,j+9))
     if j==1
         title('X (m) v/s Time(s)');
         
     end
     if j==2
         title('Y (m) v/s Time(s)');
     end
     if j==3
         title('h (m) v/s Time(s)');
         if symmetric==1
         ylim([4999,5001]);
         end
     end
         
%      if IC(j+9,1)~=0
%         axis([t0 tf 0.95*IC(j+9,1) 1.05*IC(j+9,1)]);
%      else
%         axis([t0 tf -0.1 0.1]);
%      end
     grid on
 end


figure(5);
plot(ETA)
 
figure(6)
plot3(X(:,10),X(:,11),X(:,12));
grid on
%hold on
%plot3(X(1,10),X(1,11),X(1,12),'r*');
%hold on
%plot(t,X(:,12),'r');
%hold off
% axis([0 10000 -100 100 7000 9000])

 
