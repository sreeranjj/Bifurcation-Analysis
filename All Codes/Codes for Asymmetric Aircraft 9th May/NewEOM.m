function F = RHS_func_rad_12(t0,x)

global m g Tm rho S b cbar CD0 CD_q CD_de  CL0 CL_q CL_de del_e del_a del_r eta CY_b CY_p CY_r CY_de CY_da CY_dr...
    Ixx Iyy Izz Ixz Cl_b Cl_p Cl_r Cl_de Cl_da Cl_dr Cn_b Cn_p Cn_r Cn_de Cn_da Cn_dr Cm0 Cm_q Cm_de

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



Ixy=0;
Iyz=0;
xcm=0;
ycm=0;
zcm=0;
I=eye(3);
m_o=5;
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
    

  Q=0.5*rho*(x(1)^2);                         %Dynamic Pressure     
  CD=CD0+CD_q*(cbar/2)*x(5)/x(1)+CD_de*del_e; 
  CL=CL0+CL_q*(cbar/2)*x(5)/x(1)+CL_de*del_e;
  CT=eta*Tm/Q*S;
  
 
 CX=CT-CD*cos(alpha) + CL*sin(alpha);
 CZ=-CD*sin(alpha) - CL*cos(alpha);
 CY=CY_b*x(3)+CY_p*(b/2)*x(4)/x(1)+CY_r*(b/2)*x(6)/x(1)+CY_de*del_e+CY_da*del_a+CY_dr*del_r;
 Cl=Cl_b*x(3)+Cl_p*(b/2)*x(4)/x(1)+Cl_r*(b/2)*x(6)/x(1)+Cl_de*del_e+Cl_da*del_a+Cl_dr*del_r;
 Cm=Cm0+Cm_q*(cbar/2)*x(5)/x(1)+Cm_de*del_e;
 Cn=Cn_b*x(3)+Cn_p*(b/2)*x(4)/x(1)+Cn_r*(b/2)*x(6)/x(1)+Cn_de*del_e+Cn_da*del_a+Cn_dr*del_r;
 
 
 
 mat4= [(1/m_tot)*(q*S*CX) ;
        (1/m_tot)*(q*S*CY) ;
        (1/m_tot)*(q*S*CZ) ;
             Q*S*b*Cl      ;
             Q*S*b*Cm      ;
             Q*S*b*Cn      ];
         
  
             
         
        
            
 
 
        
   
     
      
      
 




%----------------------Kinematic Equations---------------------------%
f7=x(4)+(x(5)*sin(x(7))+x(6)*cos(x(7)))*tan(x(8));
%f7=x(4)+(x(5)*sin(x(7))+x(6)*cos(x(7)))*tan(x(8));

f8=x(5)*cos(x(7))-x(6)*sin(x(7));

f9=(x(5)*sin(x(7))+x(6)*cos(x(7)))*sec(x(8));

f10=x(1)*cos(x(2))*cos(x(3))*cos(x(9))*cos(x(8))...
    +x(1)*sin(x(3))*(-sin(x(9))*cos(x(7))+sin(x(7))*sin(x(8))*cos(x(9)))...
    +x(1)*sin(x(2))*cos(x(3))*(sin(x(7))*sin(x(9))+cos(x(7))*sin(x(8))*cos(x(9)));

f11=x(1)*cos(x(2))*cos(x(3))*sin(x(9))*cos(x(8))...
    +x(1)*sin(x(3))*(cos(x(9))*cos(x(7))+sin(x(7))*sin(x(8))*sin(x(9)))...
    +x(1)*sin(x(2))*cos(x(3))*(-sin(x(7))*cos(x(9))+cos(x(7))*sin(x(8))*sin(x(9)));

f12=x(1)*(cos(x(2))*cos(x(3))*sin(x(8))-sin(x(3))*cos(x(8))*sin(x(7))...
    -sin(x(2))*cos(x(3))*cos(x(8))*cos(x(7)));