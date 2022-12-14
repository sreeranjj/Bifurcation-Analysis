function F = Equations_of_Motion(t0,State,g,Thrust,Geom,rho,cntl)

[ALPHA_BREAK,clift,clift_q,cd0,cd_q,cd_de,cy_b,croll_b,croll_p,croll_r,cm,cm_q,cn_b,cn_p,cy_r,...
    clift_de,cy_da,cy_dr,croll_da,cn_da,cn_r,cy_p,croll_dr,cn_dr,cm_de,cy_de,croll_de,cn_de] = MyAerodata;

global Velocity Altitude

% x(1)=u
% x(2)=v
% x(3)=w
% x(4)=p
% x(5)=q
% x(6)=r
% x(7)=phi
% x(8)=theta
% x(9)=psi
% x(10)=xe
% x(11)=ye
% x(12)=h

d2r = pi/180;   % deg to rad conversion
r2d = 180/pi;   % rad to deg conversion

u    =State(1);
v    =State(2);
w    =State(3);
p    =State(4);
q    =State(5);
r    =State(6);
phi  =State(7);
theta=State(8);
psi  =State(9);
xe   =State(10);
ye   =State(11);
h    =State(12);


m=Geom(1,1);
m_o=Geom(2,1);


Ixx=Geom(6,1);
Iyy=Geom(7,1);
Izz=Geom(8,1);
Ixy=Geom(9,1);
Ixz=Geom(10,1);
Iyz=Geom(11,1);

S=Geom(13,1);
b=Geom(14,1);
cbar=Geom(15,1);


xcm=Geom(18,1);
ycm=Geom(19,1);
zcm=Geom(20,1);


del_e=cntl(1);
del_a=cntl(3);
del_r=cntl(4);
eta=cntl(5);

Tm=Thrust;


 CD0=interp1(ALPHA_BREAK,cd0,alpha*r2d,'spline');
 CD_q=interp1(ALPHA_BREAK,cd_q,alpha*r2d,'spline')*r2d;
 CD_de=interp1(ALPHA_BREAK,cd_de,alpha*r2d,'spline')*r2d;
 
 CL0=interp1(ALPHA_BREAK,clift,alpha*r2d,'spline');
 CL_q=interp1(ALPHA_BREAK,clift_q,alpha*r2d,'spline')*r2d;
 CL_de=interp1(ALPHA_BREAK,clift_de,alpha*r2d,'spline')*r2d;
 
 
 
 
 CY_b=interp1(ALPHA_BREAK,cy_b,alpha*r2d,'spline')*r2d;
 CY_p=interp1(ALPHA_BREAK,cy_p,alpha*r2d,'spline')*r2d; 
 CY_r=interp1(ALPHA_BREAK,cy_r,alpha*r2d,'spline')*r2d; 
 CY_de=interp1(ALPHA_BREAK,cy_de,alpha*r2d,'spline')*r2d; 
 CY_da=interp1(ALPHA_BREAK,cy_da,alpha*r2d,'spline')*r2d; 
 CY_dr=interp1(ALPHA_BREAK,cy_dr,alpha*r2d,'spline')*r2d;
 
 
 Cl_b=interp1(ALPHA_BREAK,croll_b,alpha*r2d,'spline')*r2d; 
 Cl_p=interp1(ALPHA_BREAK,croll_p,alpha*r2d,'spline')*r2d; 
 Cl_r=interp1(ALPHA_BREAK,croll_r,alpha*r2d,'spline')*r2d; 
 Cl_de=interp1(ALPHA_BREAK,croll_de,alpha*r2d,'spline')*r2d; 
 Cl_da=interp1(ALPHA_BREAK,croll_da,alpha*r2d,'spline')*r2d; 
 Cl_dr=interp1(ALPHA_BREAK,croll_dr,alpha*r2d,'spline')*r2d;

 Cn_b=interp1(ALPHA_BREAK,cn_b,alpha*r2d,'spline')*r2d; 
 Cn_p=interp1(ALPHA_BREAK,cn_p,alpha*r2d,'spline')*r2d; 
 Cn_r=interp1(ALPHA_BREAK,cn_r,alpha*r2d,'spline')*r2d; 
 Cn_de=interp1(ALPHA_BREAK,cn_de,alpha*r2d,'spline')*r2d; 
 Cn_da=interp1(ALPHA_BREAK,cn_da,alpha*r2d,'spline')*r2d; 
 Cn_dr=interp1(ALPHA_BREAK,cn_dr,alpha*r2d,'spline')*r2d;

 Cm0=interp1(ALPHA_BREAK,cm,alpha*r2d,'spline');
 Cm_q=interp1(ALPHA_BREAK,cm_q,alpha*r2d,'spline')*r2d;
 Cm_de=interp1(ALPHA_BREAK,cm_de,alpha*r2d,'spline')*r2d;

a11=(1/m)*eta*Tm*cos(x(2))*cos(x(3));
a12=-(1/m)*0.5*rho*(x(1)^2)*S*(CD0+CD_q*(cbar/2)*x(5)/x(1)+CD_de*del_e);
a13=-g*(sin(x(8))*cos(x(2))*cos(x(3))-cos(x(8))*sin(x(7))*sin(x(3))-sin(x(2))*cos(x(3))*cos(x(7))*cos(x(8)));
f1=a11+a12+a13;

a21=x(5)-(x(4)*cos(x(2))+x(6)*sin(x(2)))*tan(x(3));
a22=-(1/m)*(1/x(1))*(1/cos(x(3)))*eta*Tm*sin(x(2));
a23=-(1/m)*(1/cos(x(3)))*0.5*rho*x(1)*S*(CL0+CL_q*(cbar/2)*x(5)/x(1)+CL_de*del_e);
a24=(1/x(1))*(1/cos(x(3)))*g*(sin(x(2))*sin(x(8))+cos(x(2))*cos(x(7))*cos(x(8)));
f2=a21+a22+a23+a24;

a31=x(4)*sin(x(2))-x(6)*cos(x(2));
a32=-(1/m)*(1/x(1))*eta*Tm*cos(x(2))*sin(x(3));
a33=-(1/m)*0.5*rho*x(1)*S*(CY_b*x(3)+CY_p*(b/2)*x(4)/x(1)+CY_r*(b/2)*x(6)/x(1)+CY_de*del_e+CY_da*del_a+CY_dr*del_r);
a34=(1/x(1))*g*(cos(x(2))*sin(x(3))*sin(x(8))+cos(x(3))*sin(x(7))*cos(x(8))-sin(x(2))*sin(x(3))*cos(x(7))*cos(x(8)));
f3=a31+a32+a33+a34;
%----------------------------------------------------------------------------------------------------------------------%

%-------------------------------------------------Moment Equations in Body Frame---------------------------------------------------------%
a41=((Iyy-Ixx-Izz)*Ixz/(Ixz^2-Ixx*Izz))*x(4)*x(5);
a42=((Ixz^2+Izz^2-Izz*Iyy)/(Ixz^2-Ixx*Izz))*x(5)*x(6);
a43=(Izz/(Ixx*Izz-Ixz^2))*0.5*rho*(x(1)^2)*S*b*(Cl_b*x(3)+Cl_p*(b/2)*x(4)/x(1)+Cl_r*(b/2)*x(6)/x(1)+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r);
a44=(Ixz/(Ixx*Izz-Ixz^2))*0.5*rho*(x(1)^2)*S*b*(Cn_b*x(3)+Cn_p*(b/2)*x(4)/x(1)+Cn_r*(b/2)*x(6)/x(1)+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r);
f4=a41+a42+a43+a44;

a51=((Izz-Ixx)/Iyy)*x(6)*x(4);
a52=-(Ixz/Iyy)*(x(4)^2-x(6)^2);
a53=(1/Iyy)*0.5*rho*(x(1)^2)*S*cbar*(Cm0+Cm_q*(cbar/2)*x(5)/x(1)+Cm_de*del_e);
f5=a51+a52+a53;

a61=((Ixx*Iyy-Ixx^2-Ixz^2)/(Ixz^2-Ixx*Izz))*x(4)*x(5);
a62=(((Ixx-Iyy+Izz)*Ixz)/(Ixz^2-Ixx*Izz))*x(5)*x(6);
a63=(Ixz/(Ixx*Izz-Ixz^2))*0.5*rho*(x(1)^2)*S*b*(Cl_b*x(3)+Cl_p*(b/2)*x(4)/x(1)+Cl_r*(b/2)*x(6)/x(1)+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r);
a64=(Ixx/(Ixx*Izz-Ixz^2))*0.5*rho*(x(1)^2)*S*b*(Cn_b*x(3)+Cn_p*(b/2)*x(4)/x(1)+Cn_r*(b/2)*x(6)/x(1)+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r);
f6=a61+a62+a63+a64;


f7=p+(q*sin(phi)+r*cos(phi))*tan(theta);
%f7=x(4)+(x(5)*sin(x(7))+x(6)*cos(x(7)))*tan(x(8));

%f8=x(5)*cos(x(7))-x(6)*sin(x(7));
f8=q*cos(phi)-r*sin(phi);


%f9=(x(5)*sin(x(7))+x(6)*cos(x(7)))*sec(x(8));
f9=(q*sin(phi)+r*cos(phi))*sec(theta);

f10=Velocity*cos(alpha)*cos(beta)*cos(psi)*cos(theta)...
    +Velocity*sin(beta)*(-sin(psi)*cos(phi)+sin(phi)*sin(theta)*cos(psi))...
    +Velocity*sin(alpha)*cos(beta)*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi));

f11=Velocity*cos(alpha)*cos(beta)*sin(psi)*cos(theta)...
    +Velocity*sin(beta)*(cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi))...
    +Velocity*sin(alpha)*cos(beta)*(-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi));

f12=Velocity*(cos(alpha)*cos(beta)*sin(theta)-sin(beta)*cos(theta)*sin(phi)...
    -sin(alpha)*cos(beta)*cos(theta)*cos(phi));


F= [f1;
    f2;
    f3;
    f4;
    f5;
    f6;
    f7;
    f8;
    f9;
    f10;
    f11;
    f12];

end