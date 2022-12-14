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

V=sqrt(u^2+v^2+w^2);
alpha=atan(w/u);
beta=asin(v/V);
I=eye(3);

m_tot= m_o + m ;

B=[  0 , zcm , -ycm ;
   -zcm ,  0  ,  xcm ;
    ycm , -xcm,   0  ];

C=[    0     , -m_tot*zcm ,  m_tot*ycm  ;
    m_tot*zcm  ,     0    ,  -m_tot*xcm ;
    -m_tot*ycm ,  m_tot*xcm ,     0      ];

D=[ Ixx , -Ixy  ,  -Ixz ;
   -Ixy ,  Iyy  ,  -Iyz ;
   -Ixz , -Iyz  ,   Izz  ];


mat=[ I , B ;
      C , D];
  
  
mat1=[                  -q*w + r*v  ;
                        -r*u + p*w  ;
                        -p*v + q*u  ;
      Ixz*p*q + (Iyy-Izz)*q*r - Ixy*r*p + Iyz*(q^2-r^2) ;
     -Iyz*p*q + Ixy*q*r + (Izz-Ixx)*r*p + Ixz*(r^2-p^2) ;
      (Ixx-Iyy)*p*q - Ixz*q*r + Iyz*r*p + Ixy*(p^2-q^2) ];
  
   
  
mat2 =[(q^2+r^2)*xcm - p*q*ycm - r*p*zcm ;
      -p*q*xcm + (r^2+p^2)*ycm - q*r*zcm ;
      -r*p*xcm - q*r*ycm + (p^2+q^2)*zcm ;
      m*(q*u-p*v)*ycm + m*(r*u-p*w)*zcm ;
      m*(p*v-q*u)*xcm + m*(r*v-q*w)*zcm ;
      m*(p*w-r*u)*xcm + m*(q*w-r*v)*ycm ];
  
  
 mat3 =[                        -g*sin(theta)       ;
                            g*cos(theta)*sin(phi)   ;
                            g*cos(theta)*cos(phi)   ;
            -zcm*m*g*cos(theta)*sin(phi) + ycm*m*g*cos(theta)*cos(phi) ;
               -zcm*m*g*sin(theta) - xcm*m*g*cos(theta)*cos(phi)        ;
                ycm*m*g*sin(theta) - xcm*m*g*cos(theta)*sin(phi)         ];
    

  Q=0.5*rho*(Velocity^2);                         %Dynamic Pressure 
  
  
  CD0=interp1(ALPHA_BREAK,cd0,alpha*r2d,'spline');
 CD_q=interp1(ALPHA_BREAK,cd_q,alpha*r2d,'spline')*r2d;
 CD_de=interp1(ALPHA_BREAK,cd_de,alpha*r2d,'spline')*r2d;
 
 CL0=interp1(ALPHA_BREAK,clift,alpha*r2d,'spline');
 CL_q=interp1(ALPHA_BREAK,clift_q,alpha*r2d,'spline')*r2d;
 CL_de=interp1(ALPHA_BREAK,clift_de,alpha*r2d,'spline')*r2d;

  
  CD=CD0+CD_q*(cbar/2)*q/Velocity+CD_de*del_e;
  D=CD*Q*S; %Drag
  

  CL=CL0+CL_q*(cbar/2)*q/Velocity+CL_de*del_e;
  L=CL*Q*S; %Lift
 
  %CT=eta*Tm/Q*S;
  T=eta*Tm;
  %Tm/(Q*S)
  %CT
  %CT*Q*S/Tm
 
  
 
 %CX=CT-CD*cos(alpha) + CL*sin(alpha);
 X=T-D*cos(alpha)+L*sin(alpha);
 %CZ=-CD*sin(alpha) - CL*cos(alpha);
 Z=-D*sin(alpha) - L*cos(alpha);
 
 %CY=CY_b*x(3)+CY_p*(b/2)*x(4)/x(1)+CY_r*(b/2)*x(6)/x(1)+CY_de*del_e+CY_da*del_a+CY_dr*del_r;
 %Cl=Cl_b*x(3)+Cl_p*(b/2)*x(4)/x(1)+Cl_r*(b/2)*x(6)/x(1)+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r;
 %Cm=Cm0+Cm_q*(cbar/2)*x(5)/x(1)+Cm_de*del_e;
 %Cn=Cn_b*x(3)+Cn_p*(b/2)*x(4)/x(1)+Cn_r*(b/2)*x(6)/x(1)+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r;
 CY_b=interp1(ALPHA_BREAK,cy_b,alpha*r2d,'spline')*r2d;
 CY_p=interp1(ALPHA_BREAK,cy_p,alpha*r2d,'spline')*r2d; 
 CY_r=interp1(ALPHA_BREAK,cy_r,alpha*r2d,'spline')*r2d; 
 CY_de=interp1(ALPHA_BREAK,cy_de,alpha*r2d,'spline')*r2d; 
 CY_da=interp1(ALPHA_BREAK,cy_da,alpha*r2d,'spline')*r2d; 
 CY_dr=interp1(ALPHA_BREAK,cy_dr,alpha*r2d,'spline')*r2d;
 
 
 Cl_b=interp1(ALPHA_BREAK,croll_b,alpha*r2d,'spline')*r2d; 
 Cl_p=interp1(ALPHA_BREAK,croll_p,alpha*r2d,'spline')*r2d; 
 Cl_r=interp1(ALPHA_BREAK,croll_r,alpha*r2d,'spline')*r2d; 
 Cl_de=interp1(ALPHA_BREAK,croll_de,alpha*r2d,'spline')*r2d ;
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
 
 


 
 CY=CY_b*beta+CY_p*(b/2)*p/Velocity+CY_r*(b/2)*r/Velocity+CY_de*del_e+CY_da*del_a+CY_dr*del_r;
 Y=CY*Q*S;
 
 Cl=Cl_b*beta+Cl_p*(b/2)*p/Velocity+Cl_r*(b/2)*r/Velocity+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r;
 clt1beta= Cl_b*beta;
 Cl_b;
 beta;
%  clt2de=Cl_de*del_e
%  Cl_de
%  del_e
%  clt3da=Cl_da*del_a
%  clt4dr=Cl_dr*del_r
 
 
% Cl=-0.007074788963297
% clt1beta =-9.717775439263039e-04
%Cl_b=-0.111903157382441
%beta=0.008684094056481
% clt2de =0
% Cl_de =0
% del_e =-0.028893927149754
% clt3da =-0.006503752066746
% clt4dr =4.007406473748893e-04
 
 
 Roll=Cl*Q*S*b;
 
 Cm=Cm0+Cm_q*(cbar/2)*q/Velocity+Cm_de*del_e;
 Cm_de*del_e;

 M=Cm*S*cbar;
 
 
 
 Cn=Cn_b*beta+Cn_p*(b/2)*p/Velocity+Cn_r*(b/2)*r/Velocity+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r;
 N=Cn*Q*S;
 
 
 mat4= [(1/m_tot)*(X) ;
        (1/m_tot)*(Y) ;
        (1/m_tot)*(Z) ;
             Q*S*b*Cl      ;
             Q*S*cbar*Cm      ;
             Q*S*b*Cn      ];
         


 BigMat=mat1+mat2+mat3+mat4;



 %MyAns=-Iyz*p*q + Ixy*q*r + (Izz-Ixx)*r*p + Ixz*(r^2-p^2) +m*(p*v-q*u)*xcm + m*(r*v-q*w)*zcm +-zcm*m*g*sin(theta) - xcm*m*g*cos(theta)*cos(phi)+Q*S*cbar*Cm 
 %ycm*m*g*sin(theta)
 %Q*S*b*Cn
 
 Fsmall=mat\BigMat;
 
 f1=Fsmall(1,1);
 f2=Fsmall(2,1);
 f3=Fsmall(3,1);
 f4=Fsmall(4,1);
 f5=Fsmall(5,1);
 f6=Fsmall(6,1);
 
 
 

%x(1)=u
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
%----------------------Kinematic Equations---------------------------%
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