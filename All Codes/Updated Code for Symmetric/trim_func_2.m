function F = trim_func_2(x)

global ALPHA_BREAK clift clift_de cd0 cd_de cm cm_de cy_de croll_de cn_de m g S rho Tm V gamma Ixx Iyy Izz Ixz b cbar

[ALPHA_BREAK,clift,clift_q,cd0,cd_q,cd_de,cy_b,croll_b,croll_p,croll_r,cm,cm_q,cn_b,cn_p,cy_r,...
    clift_de,cy_da,cy_dr,croll_da,cn_da,cn_r,cy_p,croll_dr,cn_dr,cm_de,cy_de,croll_de,cn_de] = Aerodata1;

% x(1)=V
% x(2)=alpha
% x(3)=beta
% x(4)=p
% x(5)=q
% x(6)=r
% x(7)=phi
% x(8)=theta


% x(1)=alpha
% x(2)=del_e
% x(3)=eta
%-----------------------Force Equations in Body Frame but in Terms of Wind-frame Variables--------------------------%
a11=(1/m)*x(3)*Tm*cos(x(1)*pi/180);
a12=-(1/m)*0.5*rho*(V^2)*S*(interp1(ALPHA_BREAK,cd0,x(1))+interp1(ALPHA_BREAK,cd_de,x(1))*x(2));
a13=-g*(sin((gamma+x(1))*pi/180)*cos(x(1)*pi/180)-sin(x(1)*pi/180)*cos((x(1)+gamma)*pi/180));
f1=a11+a12+a13;



a21=0;
a22=-(1/m)*(1/V)*x(3)*Tm*sin(x(1)*pi/180);
a23=-(1/m)*0.5*rho*V*S*(interp1(ALPHA_BREAK,clift,x(1))+interp1(ALPHA_BREAK,clift_de,x(1))*x(2));
a24=(1/V)*g*(sin(x(1)*pi/180)*sin((x(1)+gamma)*pi/180)+cos(x(1)*pi/180)*cos((x(1)+gamma)*pi/180));
f2=a21+a22+a23+a24;


a31=0;
a32=0;
a33=-(1/m)*0.5*rho*V*S*(interp1(ALPHA_BREAK,cy_de,x(1))*x(2));
a34=0;
f3=a31+a32+a33+a34;

a31=x(4)*sin(x(2))-x(6)*cos(x(2));
a32=-(1/m)*(1/x(1))*eta*Tm*cos(x(2))*sin(x(3));
a33=-(1/m)*0.5*rho*x(1)*S*(CY_b*x(3)+CY_p*(b/2)*x(4)/x(1)+CY_r*(b/2)*x(6)/x(1)+CY_de*del_e+CY_da*del_a+CY_dr*del_r);
a34=(1/x(1))*g*(cos(x(2))*sin(x(3))*sin(x(8))+cos(x(3))*sin(x(7))*cos(x(8))-sin(x(2))*sin(x(3))*cos(x(7))*cos(x(8)));
f3=a31+a32+a33+a34;
%----------------------------------------------------------------------------------------------------------------------%

%-------------------------------------------------Moment Equations in Body Frame---------------------------------------------------------%
a41=0;
a42=0;
a43=(Izz/(Ixx*Izz-Ixz^2))*0.5*rho*(V^2)*S*b*(interp1(ALPHA_BREAK,croll_de,x(1))*x(2));
a44=(Ixz/(Ixx*Izz-Ixz^2))*0.5*rho*(V^2)*S*b*(interp1(ALPHA_BREAK,cn_de,x(1))*x(2));
f4=a41+a42+a43+a44;

a51=0;
a52=0;
a53=(1/Iyy)*0.5*rho*(V^2)*S*cbar*(interp1(ALPHA_BREAK,cm,x(1))+interp1(ALPHA_BREAK,cm_de,x(1))*x(2));
f5=a51+a52+a53;

a61=0;
a62=0;
a63=(Ixz/(Ixx*Izz-Ixz^2))*0.5*rho*(V^2)*S*b*(interp1(ALPHA_BREAK,croll_de,x(1))*x(2));
a64=(Ixx/(Ixx*Izz-Ixz^2))*0.5*rho*(V^2)*S*b*(interp1(ALPHA_BREAK,cn_de,x(1))*x(2));
f6=a61+a62+a63+a64;
%------------------------------------------------------------------------------------------------------------------------------------------%

%----------------------Kinematic Equations---------------------------%



%--------------------------------------------------------------------%

F= [f1;
    f2;
    f3;
  ];

end
