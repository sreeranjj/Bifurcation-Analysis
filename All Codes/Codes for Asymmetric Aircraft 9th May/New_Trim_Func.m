function F = New_Trim_Func(x)

% All angles are taken in radians and trim values obtained are in radians 
 %x(1)=alpha
 %x(2)=beta
 %x(3)=eta
 %x(4)=del_e
 %x(5)=del_a
 %x(6)=del_r
    
global ALPHA_BREAK  m m_o g S rho Tm V cd0 cd_de clift clift_de cm cm_de xcm ycm zcm cy_b cy_de cy_da...
    cy_dr croll_b croll_de croll_da croll_dr cn_b cn_de cn_da cn_dr Ixx Iyy Izz Ixz b cbar

d2r=pi/180; %deg to rad conversion
r2d=180/pi; %rad to deg conversion

Ixy=0;
Iyz=0;
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
  
  

  mat3 =[                        -g*sin(x(1))      ;
                                     0             ;
                                 g*cos(x(1))       ;
                            ycm*m_tot*g*cos(x(1))  ;
               -zcm*m_tot*g*sin(x(1)) - xcm*m_tot*g*cos(x(1))        ;
                            ycm*m_tot*g*sin(x(1))         ];
                        
 
%u= V*cos(alpha)*cos(beta);
%v= V*sin(beta);
%w= V*sin(alpha)*cos(beta);

  rho
  Q=0.5*rho*(V^2)  
  %Dynamic Pressure     
  %CD=CD0 + CD_de*x(4); 
  CD= interp1(ALPHA_BREAK,cd0,x(1)*r2d,'spline') + interp1(ALPHA_BREAK,cd_de,x(1)*r2d,'spline')*x(4);
  %CL=CL0 +CL_de*x(4);
  CL=interp1(ALPHA_BREAK,clift,x(1)*r2d,'spline') + interp1(ALPHA_BREAK,clift_de,x(1)*r2d,'spline')*x(4);
  CT=x(3)*Tm/Q*S;
  
 %We are assuming flight path angle to be zero
 CX=CT-CD*cos(x(1)) + CL*sin(x(1));
 CZ=-CD*sin(x(1)) - CL*cos(x(1));
 
 %CX=CT-CD*cos(0) + CL*sin(0);
 %CZ=-CD*sin(0) - CL*cos(0);
 
 
 %CY=CY_b*x(2)+CY_de*x(4)+CY_da*x(5)+CY_dr*x(6);
 CY=interp1(ALPHA_BREAK,cy_b,x(1),'spline')*x(2)+interp1(ALPHA_BREAK,cy_de,x(1),'spline')*x(4)+interp1(ALPHA_BREAK,cy_da,x(1),'spline')*x(5)+interp1(ALPHA_BREAK,cy_dr,x(1),'spline')*x(6);
 %Cl=Cl_b*x(2)+Cl_de*x(4)+Cl_da*x(5)+Cl_dr*x(6);
 Cl=interp1(ALPHA_BREAK,croll_b,x(1),'spline')*x(2)+interp1(ALPHA_BREAK,croll_de,x(1),'spline')*x(4)+interp1(ALPHA_BREAK,croll_da,x(1),'spline')*x(5)+interp1(ALPHA_BREAK,croll_dr,x(1),'spline')*x(6);
 %Cm=Cm0+Cm_de*x(4);
 Cm=interp1(ALPHA_BREAK,cm,x(1),'spline')+interp1(ALPHA_BREAK,cm_de,x(1),'spline')*x(4);
 %Cn=Cn_b*x(2)+Cn_de*x(4)+Cn_da*x(5)+Cn_dr*x(6);
 Cn=interp1(ALPHA_BREAK,cn_b,x(1),'spline')*x(2)+interp1(ALPHA_BREAK,cn_de,x(1),'spline')*x(4)+interp1(ALPHA_BREAK,cn_da,x(1),'spline')*x(5)+interp1(ALPHA_BREAK,cn_dr,x(1),'spline')*x(6);

 
 mat4= [(1/m_tot)*(Q*S*CX) ;
        (1/m_tot)*(Q*S*CY) ;
        (1/m_tot)*(Q*S*CZ) ;
             Q*S*b*Cl      ;
             Q*S*cbar*Cm     ;
             Q*S*b*Cn      ];
         


 BigMat=mat3+mat4;
 
 F=mat\BigMat;
 
% F=[CX;
%    CY;
%    CZ;
%    Cl;
%    Cm;
%    Cn];
 
 %Here F_dot is u_dot,v_dot,w_dot,p_dot,q_dot,r_dot
 %We need alpha beta gamma
end