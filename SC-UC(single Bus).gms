* security - constraint unit commitment single Bus

Sets
t           Study Time Intervals /T0*T24/
j           number of buses /1*6/
i(j)        Generation Units /1,2,6/
vlim        voltage limit of buses /vmax,vmin/
genchar     Operation Chracteristics of Generator Units /Busno,A,B,C,Pmax,Pmin,Qmax,Qmin,Inist,MDT,MUT,Ramp,Startup,Fcost/
cc(genchar) Generation Cost Coefficients /A,B,C/
pqd         specefication for load /pd,qd/
initval     Initaial Conditions /inStat,inOnTime, inOffTime/
;

Table gendata(i,genchar) Operation Chracteristics of Units
       Busno      A      B      C      Pmax      Pmin      Qmax      Qmin      Inist      MDT      MUT      Ramp      Startup    Fcost
*                                      (MW)      (MW)     (MVAR)   (MVAR)       (h)       (h)      (h)     (MW/h)     (MBtu)    ($/MBtu)
1      1          176.9  13.5   0.0014 220       100       50        -40       4          4        4        55        100        1.2469
2      2          129.9  32.6   0.01   100       10        50        -40       2          2        2        50        200        1.2461
6      6          137.4  17.6   0.01   20        10        50        -40       2          2        2        20        0          1.2462

Parameters
pexp(cc) Exponentiation of Generated Power With repect to Cost Coefficient
**F(P(i,t))=A + B*P(i,t)**(cc) + C*P(i,t)**(cc) = F(P(i,t)) = A + B*P(i,t) + C*P(i,t)*P(i,t)
/
A 0
B 1
C 2
/
Pload(t) active Load at Each Time Interval
/
T1  219.19
T2  235.35
T3  234.67
T4  236.73
T5  239.06
T6  244.48
T7  273.39
T8  290.60
T9  283.56
T10 281.20
T11 328.61
T12 328.10
T13 326.18
T14 323.60
T15 326.86
T16 287.79
T17 260.00
T18 246.74
T19 255.97
T20 237.35
T21 243.31
T22 283.67
T23 283.05
T24 248.75
/
Qload(t) reactive Load at Each Time Interval
/
T1 50.37
T2 47.48
T3 45.62
T4 44.49
T5 44.58
T6 46.14
T7 49.85
T8 51.06
T9 53.71
T10 59.50
T11 65.73
T12 67.88
T13 69.63
T14 70.03
T15 71.55
T16 73.54
T17 73.60
T18 70.94
T19 70.72
T20 68.24
T21 68.23
T22 66.89
T23 56.33
T24 56.23
/
;
Table voltage(j,vlim)
       vmax   vmin
1      1.05   0.95
2      1.15   0.85
3      1.15   0.85
4      1.05   0.91
5      1.15   0.85
6      1.15   0.85
;
table  wind(t,j) wind power
         1 2 3 4 5     6
T1       0 0 0 0 47    0
T2       0 0 0 0 70.2  0
T3       0 0 0 0 76    0
T4       0 0 0 0 82    0
T5       0 0 0 0 84    0
T6       0 0 0 0 84    0
T7       0 0 0 0 100   0
T8       0 0 0 0 100   0
T9       0 0 0 0 78    0
T10      0 0 0 0 64    0
T11      0 0 0 0 100   0
T12      0 0 0 0 92    0
T13      0 0 0 0 84    0
T14      0 0 0 0 80    0
T15      0 0 0 0 78    0
T16      0 0 0 0 32    0
T17      0 0 0 0 4     0
T18      0 0 0 0 9     0
T19      0 0 0 0 10    0
T20      0 0 0 0 5     0
T21      0 0 0 0 6     0
T22      0 0 0 0 56    0
T23      0 0 0 0 82    0
T24      0 0 0 0 52    0
;

Variables
z          Cost Function
P(i,t)     Generated Power of i-th Unit in t-th Time Interval
Q(i,t)     Generated reactive Power of i-th Unit in t-th Time Interval
U(i,t)     On-off Status of i-th Unit in t-th Time Interval
Xon(i,t)   The Amount of Total Intervals that i-th Unit Has Been "On"  Up to t-th Interval
Xoff(i,t)  The Amount of Total Intervals that i-th Unit Has Been "Off" Up to t-th Interval
;
Positive Variable P,Xon,Xoff;
Binary Variable U;
Table Initial(i,initval) Initial state

     inStat inOnTime inOffTime
1    1      4        0
2    1      2        0
6    1      2        0
;

U.fx(i,'T0')=Initial(i,'inStat');
*Xon.fx(i,'T0')=Initial(i,'inOnTime');
*Xoff.fx(i,'T0')=Initial(i,'inOffTime');


Equations

CostFunction      Total Cost Function of Units
PLoadGen(t)       Load-Generation Balance
QLoadGen(t)       Load-Generation Balance
PUpperLimit(i,t)  Upper Limit for Power Generation
PLowerLimit(i,t)  Lower Limit for Power Generation
QUpperLimit(i,t)  Upper Limit for reactive Power Generation
QLowerLimit(i,t)  Lower Limit for reactive Power Generation
ramplimit(i,t)    Ramp limit Equation
MinUT(i,t)        Minimun Up-time   Constraint
MinDT(i,t)        Minimun Down-time Constraint
ON (i,t)          Xon  Value
OFF(i,t)          Xoff Value
;

CostFunction..                  z=e=sum(t$(ord(t)<>1),sum((i,cc),U(i,t)*gendata(i,'Fcost')*gendata(i,cc)*power(P(i,t),pexp(cc))))+sum((i,t)$(ord(t)<>1),(1-U(i,t-1))*U(i,t)*(gendata(i,'Startup')*gendata(i,'Fcost'))) ;
PLoadGen(t)$(ord(t)<>1)..       Pload(t)=e=sum(i,P(i,t)*U(i,t))+wind(t,'5');
QLoadGen(t)$(ord(t)<>1)..       Qload(t)=e=sum(i,Q(i,t)*U(i,t));
ramplimit(i,t)$(ord(t)<>1)..    P(i,t)-P(i,t-1)=l=gendata(i,'Ramp')*(1-U(i,t)*(1-U(i,t-1)))+gendata(i,'Pmin')*(U(i,t)*(1-U(i,t-1)));
PUpperLimit(i,t)$(ord(t)<>1)..  P(i,t)=l=gendata(i,'Pmax')*U(i,t);
PLowerLimit(i,t)$(ord(t)<>1)..  P(i,t)=g=gendata(i,'Pmin')*U(i,t);
QUpperLimit(i,t)$(ord(t)<>1)..  Q(i,t)=l=gendata(i,'Qmax')*U(i,t);
QLowerLimit(i,t)$(ord(t)<>1)..  Q(i,t)=g=gendata(i,'Qmin')*U(i,t);
ON (i,t)$(ord(t)<>1)..          Xon(i,t)=e=(1-U(i,t-1))*U(i,t)+U(i,t)*U(i,t-1)*(1+Xon(i,t-1));
OFF(i,t)$(ord(t)<>1)..          Xoff(i,t)=e=U(i,t-1)*(1-U(i,t))+(1-U(i,t))*(1-U(i,t-1))*(1+Xoff(i,t-1));
MinUT(i,t)$(ord(t)<>1)..       (Xon(i,t-1)-gendata(i,'MUT'))*(U(i,t-1)-U(i,t))=g=0;
MinDT(i,t)$(ord(t)<>1)..       (Xoff(i,t-1)-gendata(i,'MDT'))*(U(i,t)-U(i,t-1))=g=0;

Model  UC Unit Commitment /all/;
Option MINLP=couenne;
Solve  UC using MINLP minimizing z;
Display z.l, U.l, P.l,Q.l,Xon.l,Xoff.l, PLoadGen.m;

execute_unload "results.gdx" P.l

*=== Now write to variable levels to Excel file from GDX
*=== Since we do not specify a sheet, data is placed in first sheet
*execute 'gdxxrw.exe results.gdx o=results.xls var=U.L '

*=== Write marginals to a different sheet with a specific range
*execute 'gdxxrw.exe  results.gdx o=results.xls var=z.l rng=NewSheet'
*== Export to Excel using GDX utilities
execute 'gdxxrw.exe  results.gdx o=results.xls var=P.l'
*== Export to Excel using GDX utilities


