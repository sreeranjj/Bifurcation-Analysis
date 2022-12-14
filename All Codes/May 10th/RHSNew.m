function F = RHSNew(t0,x)

%global m g Tm rho S b cbar CD0 CD_q CD_de  CL0 CL_q CL_de del_e del_a del_r eta CY_b CY_p CY_r CY_de CY_da CY_dr...
    %Ixx Iyy Izz Ixz Cl_b Cl_p Cl_r Cl_de Cl_da Cl_dr Cn_b Cn_p Cn_r Cn_de Cn_da Cn_dr Cm0 Cm_q Cm_de

global Velocity Altitude
[F] = Equations_of_Motion(State,g,ALPHA_BREAK,F18_Aerodata,Thrust,Geom,rho,cntl);
% x(1)=V
% x(2)=alpha
% x(3)=beta
% x(4)=p
% x(5)=q
% x(6)=r
% x(7)=phi
% x(8)=theta
% x(9)=psi
% x(10)=xe
% x(11)=ye
% x(12)=h

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


m=Geom(1);
m_o=Geom(2);


Ixx=Geom(6);
Iyy=Geom(7);
Izz=Geom(8);
Ixy=Geom(9);
Ixz=Geom(10);
Iyz=Geom(11);

S=Geom(13);
b=Geom(14);
cbar=Geom(15);


xcm=Geom(18);
ycm=Geom(19);
zcm=Geom(20);


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
   -Ixz , -Ixy  ,   Izz  ];


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
      m_tot*(q*u-p*v)*ycm + m_tot*(r*u-p*w)*zcm ;
      m_tot*(p*v-q*u)*xcm + m_tot*(r*v-q*w)*zcm ;
      m_tot*(p*w-r*u)*xcm + m_tot*(q*w-r*v)*ycm ];
  
  
 mat3 =[                        -g*sin(theta)       ;
                            g*cos(theta)*sin(phi)   ;
                            g*cos(theta)*cos(phi)   ;
            -zcm*m_tot*g*cos(theta)*sin(phi) + ycm*m_tot*g*cos(theta)*cos(phi) ;
               -zcm*m_tot*g*sin(theta) - xcm*m_tot*g*cos(theta)*cos(phi)        ;
                ycm*m_tot*g*sin(theta) - xcm*m_tot*g*cos(theta)*sin(phi)         ];
    

  Q=0.5*rho*(Velocity^2);                         %Dynamic Pressure 
  
  
  CD0=interp1(ALPHA_BREAK,cd0,alpha,'spline');
 CD_q=interp1(ALPHA_BREAK,cd_q,alpha,'spline');
 CD_de=interp1(ALPHA_BREAK,cd_de,alpha,'spline');
 
 CL0=interp1(ALPHA_BREAK,clift,alpha,'spline');
 CL_q=interp1(ALPHA_BREAK,clift_q,alpha,'spline');
 CL_de=interp1(ALPHA_BREAK,clift_de,alpha,'spline');

  
  CD=CD0+CD_q*(cbar/2)*q/Velocity+CD_de*del_e;
  CL=CL0+CL_q*(cbar/2)*q/Velocity+CL_de*del_e;
  CT=eta*Tm/Q*S;
  
 
 CX=CT-CD*cos(alpha) + CL*sin(alpha);
 CZ=-CD*sin(alpha) - CL*cos(alpha);
 
 %CY=CY_b*x(3)+CY_p*(b/2)*x(4)/x(1)+CY_r*(b/2)*x(6)/x(1)+CY_de*del_e+CY_da*del_a+CY_dr*del_r;
 %Cl=Cl_b*x(3)+Cl_p*(b/2)*x(4)/x(1)+Cl_r*(b/2)*x(6)/x(1)+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r;
 %Cm=Cm0+Cm_q*(cbar/2)*x(5)/x(1)+Cm_de*del_e;
 %Cn=Cn_b*x(3)+Cn_p*(b/2)*x(4)/x(1)+Cn_r*(b/2)*x(6)/x(1)+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r;
 CY_b=interp1(ALPHA_BREAK,cy_b,alpha,'spline');
 CY_p=interp1(ALPHA_BREAK,cy_p,alpha,'spline'); 
 CY_r=interp1(ALPHA_BREAK,cy_r,alpha,'spline'); 
 CY_de=interp1(ALPHA_BREAK,cy_de,alpha,'spline'); 
 CY_da=interp1(ALPHA_BREAK,cy_da,alpha,'spline'); 
 CY_dr=interp1(ALPHA_BREAK,cy_dr,alpha,'spline');
 
 
 Cl_b=interp1(ALPHA_BREAK,croll_b,alpha,'spline'); 
 Cl_p=interp1(ALPHA_BREAK,croll_p,alpha,'spline'); 
 Cl_r=interp1(ALPHA_BREAK,croll_r,alpha,'spline'); 
 Cl_de=interp1(ALPHA_BREAK,croll_de,alpha,'spline'); 
 Cl_da=interp1(ALPHA_BREAK,croll_da,alpha,'spline'); 
 Cl_dr=interp1(ALPHA_BREAK,croll_dr,alpha,'spline');

 Cn_b=interp1(ALPHA_BREAK,cn_b,alpha,'spline'); 
 Cn_p=interp1(ALPHA_BREAK,cn_p,alpha,'spline'); 
 Cn_r=interp1(ALPHA_BREAK,cn_r,alpha,'spline'); 
 Cn_de=interp1(ALPHA_BREAK,cn_de,alpha,'spline'); 
 Cn_da=interp1(ALPHA_BREAK,cn_da,alpha,'spline'); 
 Cn_dr=interp1(ALPHA_BREAK,cn_dr,alpha,'spline');

 Cm0=interp1(ALPHA_BREAK,cm,alpha,'spline');
 Cm_q=interp1(ALPHA_BREAK,cm_q,alpha,'spline');
 Cm_de=interp1(ALPHA_BREAK,cm_de,alpha,'spline');

 
 CY=CY_b*beta+CY_p*(b/2)*p/Velocity+CY_r*(b/2)*r/Velocity+CY_de*del_e+CY_da*del_a+CY_dr*del_r;
 
 Cl=Cl_b*beta+Cl_p*(b/2)*p/Velocity+Cl_r*(b/2)*r/Velocity+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r;

 Cm=Cm0+Cm_q*(cbar/2)*q/Velocity+Cm_de*del_e;
 
 Cn=Cn_b*beta+Cn_p*(b/2)*p/Velocity+Cn_r*(b/2)*r/Velocity+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r;
 
 
 
 mat4= [(1/m_tot)*(Q*S*CX) ;
        (1/m_tot)*(Q*S*CY) ;
        (1/m_tot)*(Q*S*CZ) ;
             Q*S*b*Cl      ;
             Q*S*b*Cm      ;
             Q*S*b*Cn      ];
         


 BigMat=mat1+mat2+mat3+mat4;
 
 Fsmall=mat\BigMat;
 
 f1=Fsmall(1,1);
 f2=Fsmall(2,1);
 f3=Fsmall(3,1);
 f4=Fsmall(4,1);
 f5=Fsmall(5,1);
 f6=Fsmall(6,1);
 
 
 

%x(1)=V
% x(2)=alpha
% x(3)=beta
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