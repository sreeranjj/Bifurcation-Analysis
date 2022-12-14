function [ALPHA_BREAK,clift,clift_q,cd0,cd_q,cd_de,cy_b,croll_b,croll_p,croll_r,cm,cm_q,cn_b,cn_p,cy_r,...
    clift_de,cy_da,cy_dr,croll_da,cn_da,cn_r,cy_p,croll_dr,cn_dr,cm_de,cy_de,croll_de,cn_de] = Aerodata1


% F18 HARV aerodata (27 data points) alpha in deg
ALPHA_BREAK = [-14.0000  -10.00000   -6.00000   -2.00000    2.00000    6.00000 ...
   10.00000    14.0000    18.0000    22.0000    26.0000    30.0000    34.0000 ...
    38.0000    42.0000    46.0000    50.0000    54.0000    58.0000    62.0000 ...
    66.0000    70.0000    74.0000    78.0000    82.0000    86.0000    90.0000]';  

% All derivatives are expressed in /deg

% here is lift-coefficient due to alfa
clift = [-1.03860  -0.880907  -0.551318  -0.216206   0.154461   0.560857 ...
	     0.951575    1.24168    1.39786    1.53103    1.64110    1.77724    1.88652 ...
	     1.88940    1.80975    1.68556    1.57524    1.46011    1.31884    1.17117 ...
	     1.02151   0.858875   0.685334   0.529482   0.378268   0.229349    9.10933E-02];

% here is clift_q (lift coeff due to pitch rate)
clift_q = [7.87143E-02    7.87143E-02    7.87144E-02    7.88452E-02    7.64891E-02 ...
	     6.91587E-02    6.33554E-02    5.99084E-02    5.89049E-02    6.05629E-02 ...
	     6.28755E-02    6.91587E-02    8.89245E-02    0.119991   0.147829   0.159698 ...
	     0.131772       1.03847E-01    9.21534E-02    8.79646E-02    8.63938E-02    8.63938E-02 ...
	     8.49975E-02    8.25541E-02    7.97615E-02    7.66200E-02    7.24312E-02]; 

% here is cd0 (zero-lift drag coeff)
cd0 =  [0.227325   0.170673    7.47685E-02    3.20977E-02    2.79305E-02 ...
	     6.13288E-02   0.154023   0.261263   0.363863   0.486635   0.643989   0.863792 ...
	     1.08558    1.27900    1.47638    1.61911    1.74604    1.84448    1.92402 ...
	     1.98384    2.02834    2.06708    2.11362    2.15514    2.17363    2.16283 ...
	     2.12822];

% here is cd_q (drag coeff due to pitch rate)
cd_q = [ 1.48740E-16   -2.21604E-16   -8.70391E-05   -7.79330E-04   -8.31959E-04 ...
	    -6.60720E-05    1.33279E-03    1.67685E-04    5.84481E-06   -9.59697E-17 ...
	    -6.48262E-17   -2.05428E-16   -2.53047E-16    2.02363E-16    3.51664E-16 ...
	     8.52635E-16    5.44982E-16   -5.45460E-16    1.14528E-15   -2.50755E-16 ...
	     8.43503E-16   -8.45493E-16   -2.54842E-15   -7.66583E-16   -3.79368E-16 ...
	     9.11688E-16    9.07194E-16];

% here is cd_del (drag coeff due to left stabilator deflection)
cd_del = [ -2.63851E-03   -2.42289E-03   -2.00960E-03   -1.61668E-03   -1.11178E-03 ...
	    -4.31951E-04    4.33892E-04    1.11520E-03    1.97203E-03    2.82682E-03 ...
	     3.59917E-03    4.29851E-03    4.93824E-03    5.55955E-03    5.88711E-03 ...
	     5.82388E-03    5.49440E-03    5.22574E-03    5.15071E-03    5.18088E-03 ...
	     5.20101E-03    5.06911E-03    3.65306E-03    2.50369E-03    1.69598E-03 ...
	     1.25007E-03    1.13131E-03];

% here is cd_der (drag coeff due to right stabilator deflection)
cd_der = [ -2.63851E-03   -2.42289E-03   -2.00960E-03   -1.61668E-03   -1.11178E-03 ...
	    -4.31951E-04    4.33892E-04    1.11520E-03    1.97203E-03    2.82682E-03 ...
	     3.59917E-03    4.29851E-03    4.93824E-03    5.55955E-03    5.88711E-03 ...
	     5.82388E-03    5.49440E-03    5.22574E-03    5.15071E-03    5.18088E-03 ...
	     5.20101E-03    5.06911E-03    3.65306E-03    2.50369E-03    1.69598E-03 ...
	     1.25007E-03    1.13131E-03];

cd_de = cd_del+cd_der;

% here is sideforce coeff due to sideslip derivative
cy_b = [-1.71957E-02   -1.71957E-02   -1.77407E-02   -1.82857E-02   -1.84197E-02 ...
	   -1.84159E-02   -1.87763E-02   -1.77142E-02   -1.55057E-02   -1.40980E-02 ...	
	   -1.28732E-02   -1.26852E-02   -1.25160E-02   -1.35594E-02   -1.40314E-02 ...
	   -1.33128E-02   -1.19496E-02   -1.30678E-02   -1.37112E-02   -1.37528E-02 ...
	   -1.33426E-02   -1.29076E-02   -1.33076E-02   -1.38320E-02   -1.37811E-02 ...
	   -1.31714E-02   -1.27360E-02];

% here is croll_b (rolling moment coeff due to sideslip)
croll_b = [-3.70196E-05   -3.70196E-05   -4.41520E-04   -8.46019E-04   -1.27065E-03 ...
	   -1.70123E-03   -2.11290E-03   -2.77077E-03   -3.37291E-03   -3.47757E-03 ...
	   -3.16486E-03   -2.35754E-03   -1.22554E-03   -1.14567E-03   -1.63756E-03 ...
	   -2.22957E-03   -2.45882E-03   -2.47686E-03   -2.44493E-03   -2.51646E-03 ...
	   -2.70976E-03   -2.90803E-03   -2.82658E-03   -2.75756E-03   -2.70552E-03 ...
	   -2.67292E-03   -2.66024E-03]; 

% here is croll_p (rolling moment coeff due to roll rate)
croll_p = [-7.05113E-03   -7.05113E-03   -7.05113E-03   -7.05113E-03   -7.05113E-03 ...
	   -7.05113E-03   -6.24828E-03   -5.25344E-03   -4.62512E-03   -3.85718E-03 ...
	   -5.38434E-03   -7.94125E-03   -8.68301E-03   -8.72665E-03   -4.88692E-03 ...
	    3.49064E-05   -3.31613E-03   -3.73500E-03   -4.04916E-03   -4.36332E-03 ...
	   -4.67748E-03   -4.88692E-03   -5.02655E-03   -5.16617E-03   -5.30580E-03 ...
	   -5.41052E-03   -5.41052E-03];

% here is croll_r (rolling moment coeff due to yaw rate)
croll_r = [2.14676E-04    2.14676E-04    2.14676E-04    5.85558E-04    1.52804E-03 ...
	    2.57087E-03    3.39117E-03    3.84496E-03    3.98895E-03    4.11549E-03 ...
	    4.43401E-03    5.18886E-03    5.68192E-03    5.06669E-03    4.07360E-03 ...
	    3.28122E-03    2.56214E-03    2.18515E-03    1.86052E-03    1.58127E-03 ...
	    1.33692E-03    1.11352E-03    1.12748E-03    1.13097E-03    9.28515E-04 ...
	    5.34071E-04    1.71042E-04];

% here is cm (pitching moment coeff due to angle-of-attack)
cm = [9.62520E-02    8.07290E-02    5.54066E-02    2.80686E-02    5.08503E-03 ...
	    -6.61034E-03   -1.64371E-02   -1.60245E-02   -5.86329E-02   -8.42348E-02 ...
	    -8.52663E-02  -0.105573  -0.117061  -0.110115  -0.114069  -0.118398   -1.01818E-01 ...
	   -0.146057  -0.228343  -0.351592  -0.426319  -0.476887  -0.601749  -0.563499 ...
	   -0.547688  -0.570013  -0.588038]; 

% here is cm_q (pitching moment coeff due to pitch rate)
cm_q = [-8.68301E-02   -8.68301E-02   -8.68301E-02   -8.65683E-02   -8.53902E-02 ...
	    -8.08524E-02   -7.57473E-02   -7.23875E-02   -7.56600E-02   -8.29032E-02 ...
	    -8.81828E-02   -1.02364E-01  -0.147655  -0.205556  -0.132732    1.80642E-02 ...
	     4.31969E-02    1.80642E-02   -1.96350E-02   -6.50135E-02   -1.00967E-01 ...
	    -9.81748E-02   -9.53822E-02   -9.25897E-02   -8.90991E-02   -8.52593E-02 ...
	    -8.24668E-02];

% here is cn_b (yawing moment coeff due to sideslip)
cn_b = [1.42071E-03    1.42071E-03    1.52521E-03    1.62971E-03    1.65089E-03 ...
	    1.67084E-03    1.80024E-03    1.46016E-03    9.91781E-04    6.82316E-04 ...
	    1.75902E-04   -1.09086E-03   -1.11070E-03   -1.59428E-03   -1.88776E-03 ...
	   -1.76721E-03   -1.35917E-03   -1.52852E-03   -1.63273E-03   -1.59147E-03 ...
	   -1.44456E-03   -1.34440E-03   -1.47796E-03   -1.60891E-03   -1.69618E-03 ...
	    -1.72843E-03   -1.73807E-03];

% here is cn_p (yawing moment coeff due to roll rate)
cn_p = [-1.26406E-03   -1.26406E-03   -1.26406E-03   -1.23533E-03   -1.17787E-03 ...
	    -1.12041E-03   -9.39611E-04   -6.36739E-04   -2.47236E-04    1.00356E-04 ...
	     3.70882E-04    4.58149E-04    3.97062E-04    2.48709E-04    2.33874E-04 ...
	     4.18879E-04    6.98132E-04    9.77384E-04   -1.04720E-03   -1.60570E-03 ...
	    -3.24631E-04   -2.26893E-04   -1.44164E-03   -6.98132E-04   -8.81388E-14 ...
	    -9.39512E-14   -9.97960E-14]; 

% here is cn_r (yawing moment coeff due to yaw rate)
cn_r = [-3.11978E-03   -3.11978E-03   -3.11978E-03   -3.11105E-03   -3.10669E-03 ...
	    -3.12850E-03   -3.17650E-03   -3.27686E-03   -3.46448E-03   -3.97935E-03 ...
	    -4.62949E-03   -5.14872E-03   -5.53269E-03   -5.63741E-03   -5.53706E-03 ...
	    -5.23162E-03   -4.53349E-03   -3.97499E-03   -2.89288E-03   -1.70606E-03 ...
	    -6.93768E-04    4.36324E-06   -6.93769E-04   -2.64854E-03   -4.18443E-03 ...
	    -4.81274E-03   -5.23162E-03];

% here is cy_p (sideforce coeff due to roll rate)
cy_p = [-1.06465E-03   -1.06465E-03   -1.06465E-03   -8.63938E-04   -2.09440E-04 ...
	    4.14516E-04    5.80322E-04    5.32325E-04    4.53786E-04    2.22529E-04 ...
	   -1.74533E-04   -5.49779E-04   -9.42478E-04   -1.35699E-03   -2.60752E-03 ...
	   -4.39823E-03   -5.23599E-03   -3.00197E-03   -1.65806E-03   -7.85399E-04 ...
	   -2.96706E-04   -4.36332E-04   -9.25024E-04   -9.94838E-04   -9.25024E-04 ...
	   -9.77384E-04   -1.39626E-03];

% here is cy_r (sideforce coeff due to yaw rate)
cy_r = [2.69653E-03    2.69653E-03    2.69653E-03    2.89725E-03    3.36412E-03 ...
	    3.94444E-03    4.45932E-03    3.82227E-03    1.22609E-03   -1.61443E-03 ...
	   -3.46448E-03   -4.66876E-03   -5.36689E-03   -5.67232E-03   -3.89208E-03 ...
	    3.49069E-05    5.06146E-03    5.06145E-03    4.11898E-03    2.79253E-03 ...
	    1.60570E-03    1.04720E-03    8.37758E-04    6.28319E-04    4.88692E-04 ...
	    4.18879E-04    3.49066E-04];

% Control derivatives begin here 

% here is clift_del (lift coeff due to left stabilator deflection)
clift_del = [ 5.03735E-03    5.30820E-03    5.80124E-03    6.25003E-03    6.76123E-03 ...
	     7.18902E-03    7.15986E-03    6.87948E-03    6.62941E-03    6.28904E-03 ...
	     5.97515E-03    5.90878E-03    5.86386E-03    5.55871E-03    4.76025E-03 ...
	     3.56060E-03    2.26618E-03    1.55859E-03    1.04476E-03    5.88013E-04 ...
	     1.18163E-04   -3.68531E-04   -1.10325E-03   -1.79305E-03   -2.27413E-03 ...
	    -2.48706E-03   -2.47676E-03];
 
% here is clift_der (lift coeff due to right stabilator deflection)
clift_der = [ 5.03735E-03    5.30820E-03    5.80124E-03    6.25003E-03    6.76123E-03 ...
	     7.18902E-03    7.15986E-03    6.87948E-03    6.62941E-03    6.28904E-03 ...
	     5.97515E-03    5.90878E-03    5.86386E-03    5.55871E-03    4.76025E-03 ...
	     3.56060E-03    2.26618E-03    1.55859E-03    1.04476E-03    5.88013E-04 ...
	     1.18163E-04   -3.68531E-04   -1.10325E-03   -1.79305E-03   -2.27413E-03 ...
	    -2.48706E-03   -2.47676E-03];

% here is clift_de
clift_de = clift_del+clift_der;

% here is cy_da (sideforce coeff due to aileron deflection)
cy_da = [-5.77397E-04   -5.77397E-04   -5.77397E-04   -5.46240E-04   -4.60488E-04 ...
	   -4.06750E-04   -4.16183E-04   -3.70449E-04   -1.95515E-04    4.51627E-05 ...
	    2.63830E-04    4.73923E-04    6.38852E-04    7.09741E-04    6.50343E-04 ...
	    1.99345E-04   -1.00616E-04   -1.48179E-04   -1.74477E-04   -1.95743E-04 ...
	   -2.18725E-04   -2.40677E-04   -2.78179E-04   -3.06420E-04   -3.28602E-04 ...
	   -3.47581E-04   -3.65875E-04];

% here is cy_dr (sideforce coeff due to rudder deflection)
cy_dr = [3.46480E-03    3.46480E-03    3.46480E-03    3.52987E-03    3.66813E-03 ...
	    3.73320E-03    3.66813E-03    3.46480E-03    3.08253E-03    2.57827E-03 ...
	    2.11467E-03    1.79747E-03    1.59413E-03    1.47213E-03    1.50955E-03 ...
	    1.59771E-03    1.48189E-03    1.25806E-03    1.04497E-03    9.72747E-04 ...
	    1.04790E-03    1.15331E-03    1.10125E-03    1.01504E-03    1.04757E-03 ...
	    1.15266E-03    1.08499E-03];

% here is croll_da (rolling moment coeff due to ailreon deflection)
croll_da = [1.16543E-03    1.16543E-03    1.16543E-03    1.16686E-03    1.17021E-03 ...
	    1.17426E-03    1.14013E-03    1.04034E-03    8.40294E-04    6.31892E-04 ...
	    5.23513E-04    4.56671E-04    3.86010E-04    3.20124E-04    2.81594E-04 ...
	    2.56481E-04    2.33946E-04    2.11792E-04    1.97660E-04    1.86966E-04 ...
	    1.74839E-04    1.56123E-04    9.61564E-05    5.02267E-05    1.77607E-05 ...
	   -6.65660E-23   -7.15216E-23];

% here is cn_da (yawing moment coeff due to aileron deflection)
cn_da = [-8.69999E-06   -8.69999E-06   -8.69999E-06   -2.06625E-05   -3.45825E-05 ...
	    -3.93675E-05   -4.78500E-05   -6.15525E-05   -7.74300E-05   -1.02660E-04 ...
	    -1.32240E-04   -1.55512E-04   -1.72042E-04   -1.80090E-04   -1.97403E-04 ...
	    -2.05320E-04   -2.40120E-04   -2.69352E-04   -2.87796E-04   -3.02644E-04 ...
	    -1.09272E-04    2.45624E-12    8.35200E-05    1.59210E-04    2.21850E-04 ...
	     2.43600E-04    1.74000E-04];

% here is croll_dr (rolling moment coeff due to rudder deflection)
croll_dr = [2.86293E-04    2.86293E-04    2.86293E-04    2.69213E-04    2.39120E-04 ...
	    2.18787E-04    2.10653E-04    2.07400E-04    1.80560E-04    1.24440E-04 ...
	    6.91333E-05    3.49733E-05    1.87067E-05    4.88000E-06   -3.44528E-05 ...
	   -8.67339E-05   -1.08661E-04   -7.83403E-05   -2.63520E-05    3.68277E-05 ...
	    9.48021E-05    1.25253E-04    9.58432E-05    7.66811E-05    6.73115E-05 ...
	    6.47088E-05    6.32773E-05];

% here is cn_dr (yawing moment coeff due to rudder deflection)
cn_dr = [-1.15493E-03   -1.15493E-03   -1.15493E-03   -1.16632E-03   -1.18177E-03 ...
	    -1.18259E-03   -1.15656E-03   -1.08336E-03   -9.58920E-04   -8.03573E-04 ...
	    -6.80760E-04   -6.18947E-04   -5.90480E-04   -5.80720E-04   -5.15533E-04 ...
	    -4.17092E-04   -3.72355E-04   -3.86914E-04   -3.85538E-04   -3.51603E-04 ...
	    -3.00387E-04   -2.75612E-04   -2.70873E-04   -3.22200E-04   -3.84710E-04 ...
	    -4.23059E-04   -4.16915E-04];

% here is cm_del (pitching moment coeff due to left stabilator deflection)
cm_del = [-6.46703E-03   -6.58967E-03   -6.83004E-03   -7.09553E-03   -7.45696E-03 ...
	    -7.82965E-03   -8.00920E-03   -8.04445E-03   -8.05054E-03   -8.00358E-03 ...
	    -7.99808E-03   -7.98872E-03   -7.79706E-03   -7.35454E-03   -6.39032E-03 ...
	    -5.05660E-03   -3.75061E-03   -3.37707E-03   -3.03092E-03   -2.32363E-03 ...
	    -1.37503E-03   -8.13315E-04   -7.80549E-04   -1.77291E-03   -2.35148E-03 ...
	    -2.35499E-03   -2.89985E-03];
 
% here is cm_der (pitching moment coeff due to right stabilator deflection)
cm_der = [ -6.46703E-03   -6.58967E-03   -6.83004E-03   -7.09553E-03   -7.45696E-03 ...
	    -7.82965E-03   -8.00920E-03   -8.04445E-03   -8.05054E-03   -8.00358E-03 ...
	    -7.99808E-03   -7.98872E-03   -7.79706E-03   -7.35454E-03   -6.39032E-03 ...
	    -5.05660E-03   -3.75061E-03   -3.37707E-03   -3.03092E-03   -2.32363E-03 ...
	    -1.37503E-03   -8.13315E-04   -7.80549E-04   -1.77291E-03   -2.35148E-03 ...
	    -2.35499E-03   -2.89985E-03];

cm_de = cm_del+cm_der;

% here is cy_del (sideforce due to left stabilator deflection) 
cy_del =[-1.21826E-03   -1.21826E-03   -1.21826E-03   -1.18543E-03   -1.08596E-03 ...
	   -9.69152E-04   -8.62255E-04   -7.02810E-04   -4.98886E-04   -2.88800E-04 ...
	   -4.99497E-05    2.06273E-04    4.27489E-04    5.91445E-04    7.25571E-04 ...
	    8.56568E-04    1.00528E-03    1.13120E-03    1.22832E-03    1.30397E-03 ...
	    1.29672E-03    1.07641E-03    5.04952E-04   -4.54204E-05   -5.69973E-04 ...
	   -1.00132E-03   -1.20941E-03];

% here is cy_der (sideforce due to right stabilator deflection) 
cy_der = [1.21826E-03    1.21826E-03    1.21826E-03    1.18543E-03    1.08596E-03 ...
	    9.69152E-04    8.62255E-04    7.02810E-04    4.98886E-04    2.88800E-04 ...
	    4.99497E-05   -2.06273E-04   -4.27489E-04   -5.91445E-04   -7.25571E-04 ...
	   -8.56568E-04   -1.00528E-03   -1.13120E-03   -1.22832E-03   -1.30397E-03 ...
	   -1.29672E-03   -1.07641E-03   -5.04952E-04    4.54205E-05    5.69973E-04 ...
	    1.00132E-03    1.20941E-03];

cy_de = cy_del+cy_der;

% here is croll_del (rolling moment due to left stabilator deflection)
croll_del = [5.65136E-04    5.65136E-04    5.65136E-04    5.80368E-04    6.12176E-04 ...
	    6.47732E-04    6.78982E-04    6.89279E-04    6.85958E-04    6.90546E-04 ...
	    6.99129E-04    6.97942E-04    6.82421E-04    6.50934E-04    5.94802E-04 ...
	    5.18515E-04    4.32725E-04    3.44090E-04    2.67459E-04    2.01739E-04 ...
	    1.42504E-04    8.19937E-05    2.43887E-05   -5.32750E-06   -2.92725E-05 ...
	   -5.58400E-05   -7.97000E-05];

% here is croll_der (rolling moment due to right stabilator deflection)
croll_der = [-5.65136E-04   -5.65136E-04   -5.65136E-04   -5.80368E-04   -6.12176E-04 ...
	   -6.47732E-04   -6.78982E-04   -6.89279E-04   -6.85958E-04   -6.90546E-04 ...
	   -6.99129E-04   -6.97942E-04   -6.82421E-04   -6.50934E-04   -5.94802E-04 ...
	   -5.18515E-04   -4.32725E-04   -3.44090E-04   -2.67459E-04   -2.01739E-04 ...
	   -1.42504E-04   -8.19937E-05   -2.43887E-05    5.32750E-06    2.92725E-05 ...
	    5.58400E-05    7.97000E-05]; 

croll_de = croll_del+croll_der;

% here is cn_del (yawing moment due to left stabilator deflection)
cn_del = [ 4.45020E-04    4.45020E-04    4.45020E-04    4.36729E-04    4.17690E-04 ...
	     3.94303E-04    3.58695E-04    2.82299E-04    1.66782E-04    6.38125E-05 ...
	    -7.86875E-06   -9.12188E-05   -2.04675E-04   -3.33137E-04   -4.30339E-04 ...
	    -4.76974E-04   -4.53709E-04   -5.23318E-04   -6.41640E-04   -7.97428E-04 ...
	    -9.85781E-04   -1.20466E-03   -1.53810E-03   -1.97276E-03   -2.31246E-03 ...
	    -2.51679E-03   -2.70661E-03];
 
% here is cn_der (yawing moment due to right stabilator deflection)
cn_der = [ -4.45020E-04   -4.45020E-04   -4.45020E-04   -4.36729E-04   -4.17690E-04 ...
	    -3.94303E-04   -3.58695E-04   -2.82299E-04   -1.66782E-04   -6.38125E-05 ...
	     7.86875E-06    9.12188E-05    2.04675E-04    3.33137E-04    4.62771E-04 ...
	     5.89670E-04    7.57910E-04    8.80937E-04    1.03114E-03    1.17862E-03 ...
	     1.32721E-03    1.51793E-03    1.84599E-03    2.19195E-03    2.41929E-03 ...
	     2.53951E-03    2.70661E-03];

cn_de = cn_del+cn_der;

