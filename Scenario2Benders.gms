Set
i     Number of Buses    /1*6/
pv(i) pv buses  /1,2,4/
pq(i) pq buses /3,5,6/
t     Time Intervals  /t1*t24/
iter  iteration of algorithm /iter1*iter18/;
Alias(i,j);
*====================================================================================================================================================================================
Table Pload(t,i)  Exist PLoads At Buses In 24 Time Intervals
    3
t1  219.19
t2  235.35
t3  234.67
t4  236.73
t5  239.06
t6  244.48
t7  273.39
t8  290.60
t9  283.56
t10 281.20
t11 328.61
t12 328.10
t13 326.18
t14 323.60
t15 326.86
t16 287.79
t17 260.00
t18 246.74
t19 255.97
t20 237.35
t21 243.31
t22 283.67
t23 283.05
t24 248.75
;
*--------------------------------------------------------------------------------------------------------------------------------------------------------------
Table Gendata(pv,*) Generation Data

    C    pmax  pmin   ST

1   13   220   90     300
2   32   100   10     200
4   17   20    10     200 ;
*------------------------------------------------------------------------------------------------------
table  wind(t,j) wind power
         1 2 3 4       5  6
t1       0 0 0 47.9    0  0
t2       0 0 0 73.1    0  0
t3       0 0 0 73.8    0  0
t4       0 0 0 73      0  0
t5       0 0 0 81.2    0  0
t6       0 0 0 81      0  0
t7       0 0 0 105.7   0  0
t8       0 0 0 99.4    0  0
t9       0 0 0 85.1    0  0
t10      0 0 0 62.1    0  0
t11      0 0 0 101.2   0  0
t12      0 0 0 92.2    0  0
t13      0 0 0 85.8    0  0
t14      0 0 0 78.5    0  0
t15      0 0 0 81.8    0  0
t16      0 0 0 31      0  0
t17      0 0 0 4.3     0  0
t18      0 0 0 8       0  0
t19      0 0 0 9.8     0  0
t20      0 0 0 4.6     0  0
t21      0 0 0 5.6     0  0
t22      0 0 0 50.7    0  0
t23      0 0 0 74.7    0  0
t24      0 0 0 49.7    0  0
;
;
*--------------------------------------------------------------------------------------------------------------------------------------------------------------
 Table Pf_Max(i,j)   Maximum Allowed Power Folw Between Buses

      1        2        3     4
1     0        100      100   100
2     100      0        100    0
3     200      100      0     100
4     100      0        100    0
5     0        0        200    0
6     8        0        100    0
;
*=============================================================================================================================================================
scalar
BaseMVA /100/
*------------------------------------------------------------------------------------------------------------------------------------------------------------------
Binary Variable
U(i,t)                      On Or Off Status Of Unit i At Time t
S(pv,t)
s_optcut
s_fcut;
U.l(pv,t)=1;
U.l(pq,t)=0;
S.l(pv,t)=0;
s_optcut.l=0;
s_fcut.l=0;;

*/////////////////////////////////////////////////////////// Benders master problem ////////////////////////////////////////////////////////////////
sets
optcut(iter)  'dynamic optimal cut set'
fcut(iter)  'dynamic optimal cut set'
;
optcut(iter)=no;
fcut(iter)=no;

scalar LC /1/;
*-----------------------------------------------------------------------------------------------------------------------------------------------------------------
variables
z_master 'the objective function of main problem (lower bound of original problem)'
z_opt
z_opt1
;
parameter
U_(i,t,iter)
S_(pv,t,iter)
pma_(pv,t,iter)
pmi_(pv,t,iter)
vt_(t,iter)
pma_check_(pv,t,iter)
pmi_check_(pv,t,iter)
;
pma_(pv,t,iter)=0;
pmi_(pv,t,iter)=0;
pma_check_(pv,t,iter)=0;
pmi_check_(pv,t,iter)=0;
vt_(t,iter)=0;
*==============================================================================================================================================================================================================================================================================
equations
Master_obj_eq

Le_BtoB_1(pv,t)
Le_BtoB_2(pv,t)
Le_BtoB_3(pv,t)

master_reserve(t)          DC Power Flow

eq1(pq,t)
;

Master_obj_eq            ..  z_master=g=sum((pv,t),(U(pv,t)-S(pv,t))*Gendata(pv,'ST'));

Le_BtoB_1(pv,t)          ..  S(pv,t)=l=U(pv,t);
Le_BtoB_2(pv,t)          ..  S(pv,t)=l=U(pv,t-1);
Le_BtoB_3(pv,t)          ..  S(pv,t)=g=U(pv,t)+(U(pv,t-1)-1);

master_reserve(t)        ..  sum(pv,Gendata(pv,'pmax')*U(pv,t))+sum(i,wind(t,i))=g=sum(i,Pload(t,i));

eq1(pq,t)                ..  U(pq,t)=e=0;
*==============================================================================================================================================================================================================================================================================
model master_problem_ini /Master_obj_eq,Le_BtoB_1,Le_BtoB_2,Le_BtoB_3,master_reserve,eq1/;

equation
feasiblecut_eq(t,iter);


feasiblecut_eq(t,fcut)             ..  vt_(t,fcut)+sum(pv,pma_check_(pv,t,fcut)*Gendata(pv,'pmax')*(U(pv,t)-U_(pv,t,fcut))-pmi_check_(pv,t,fcut)*Gendata(pv,'pmin')*(U(pv,t)-U_(pv,t,fcut)))=l=0;

model master_problem_fcut /Master_obj_eq,feasiblecut_eq,Le_BtoB_1,Le_BtoB_2,Le_BtoB_3,master_reserve,eq1/;

equation
optimalitycut_eq(iter);


optimalitycut_eq(optcut)            ..  z_master=g=sum((pv,t),(U(pv,t)-S(pv,t))*Gendata(pv,'ST'))+z_opt1.l+sum((pv,t),pma_(pv,t,optcut)*Gendata(pv,'pmax')*(U(pv,t)-U_(pv,t,optcut))-pmi_(pv,t,optcut)*Gendata(pv,'pmin')*(U(pv,t)-U_(pv,t,optcut)));

model master_problem_optcut /Master_obj_eq,optimalitycut_eq,Le_BtoB_1,Le_BtoB_2,Le_BtoB_3,master_reserve,eq1/;

model master_problem /Master_obj_eq,optimalitycut_eq,feasiblecut_eq,Le_BtoB_1,Le_BtoB_2,Le_BtoB_3,master_reserve,eq1/;

*////////////////////////////////////////////////////////////// feasiblity_check //////////////////////////////////////////////////////////////////
Variables
z_ckeck                     Objective Variable
feed_check(i,j,t)
;
*-----------------------------------------------------------------------------------------------------------------------------------------------------------------
Positive Variables
Pg_check(i,t)                     Generation Active Power Of Unit i At Time t
LS(i,t)
LS_sum
vt(t);

Pg_check.fx(pq,t)=0;
LS.up(i,t)=Pload(t,i);
LS_sum.l=0;
*==================================================================================================================================================================
Equations
Obj_ckeck                         Objective Function That Containes Total Cost Of Generation
Power_Flow_ckeck(i,t)             DC Power Flow
Line_eq_check(i,j,t)
LS_sum_eq
vt_eq(t)
*Line_Limit1_ckeck(i,j,t)          Maximum Allowed Power Folw Between Buses Function1
*Line_Limit2_ckeck(i,j,t)          Maximum Allowed Power Folw Between Buses Function2
pmi_check(pv,t)
pma_check(pv,t)
;
*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Obj_ckeck                  ..  z_ckeck=e=sum((i,t),LS(i,t)*LC);
Power_Flow_ckeck(i,t)      ..  Pg_check(i,t)+LS(i,t)-Pload(t,i)+wind(t,i)=e=sum(j,feed_check(i,j,t));
Line_eq_check(i,j,t)       ..  feed_check(i,j,t)=e=-feed_check(j,i,t);

LS_sum_eq                  ..  LS_sum=e=sum((i,t),LS(i,t));
vt_eq(t)                   ..  vt(t)=e=sum(i,LS(i,t));

*Line_Limit1_ckeck(i,j,t)   ..  feed_check(i,j,t)=l=Pf_Max(i,j);
*Line_Limit2_ckeck(i,j,t)   ..  feed_check(i,j,t)=g=-Pf_Max(i,j);

pmi_check(pv,t)            ..   -Pg_check(pv,t)=l=-Gendata(pv,'pmin')*U.l(pv,t);
pma_check(pv,t)            ..   Pg_check(pv,t)=l=Gendata(pv,'pmax')*U.l(pv,t);

model feasiblity_check /Obj_ckeck,Power_Flow_ckeck,Line_eq_check,LS_sum_eq,vt_eq,pmi_check,pma_check/;

*/////////////////////////////////////////////////////////////// optimal_operation /////////////////////////////////////////////////////////////////

Variables
feed(i,j,t)
;
*-----------------------------------------------------------------------------------------------------------------------------------------------------------------
Positive Variables
Pg(i,t)                     Generation Active Power Of Unit i At Time t
;

Pg.fx(pq,t)=0;
*==================================================================================================================================================================
Equations
obj_eq
Obj                         Objective Function That Containes Total Cost Of Generation
Power_Flow(i,t)             DC Power Flow
Line_eq(i,j,t)
*Line_Limit1(i,j,t)          Maximum Allowed Power Folw Between Buses Function1
*Line_Limit2(i,j,t)          Maximum Allowed Power Folw Between Buses Function2
pmi(pv,t)
pma(pv,t)
;
*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Obj                      ..  z_opt=e=z_opt1+sum((pv,t),(U.l(pv,t)-S.l(pv,t))*Gendata(pv,'ST'));
obj_eq                   ..  z_opt1=e=sum((pv,t),Pg(pv,t)*Gendata(pv,'C'));
Power_Flow(i,t)          ..  Pg(i,t)-Pload(t,i)+wind(t,i)=e=sum(j,feed(i,j,t));

Line_eq(i,j,t)           ..  feed(i,j,t)=e=-feed(j,i,t);

*Line_Limit1(i,j,t)       ..  feed(i,j,t)=l=Pf_Max(i,j);
*Line_Limit2(i,j,t)       ..  feed(i,j,t)=g=-Pf_Max(i,j);

pmi(pv,t)                ..   -Pg(pv,t)=l=-Gendata(pv,'pmin')*U.l(pv,t);
pma(pv,t)                ..   Pg(pv,t)=l=Gendata(pv,'pmax')*U.l(pv,t);

model optimal_operation /Obj,obj_eq,Power_Flow,Line_eq,pmi,pma/;
*///////////////////////////////////////////////////////////////// benders algorithm ///////////////////////////////////////////////////////////////
parameter results(iter,*);
parameter nonconverged /yes/;
scalar  unbounded  /1.0e9/;
scalar uperbound /+inf/;
scalar lowerbound /-inf/;

option LP=CPLEX;
option threads=4;
solve  optimal_operation minimizing z_opt using lp;
uperbound=z_opt.l;



loop(iter$(nonconverged),

option MIP=CPLEX;
option optcr=0.01;
option threads=4;

if (s_optcut.l=0 and s_fcut.l=0,
solve master_problem_ini minimizing z_master using MIP;
display z_master.l,U.l, S.l, Pg.l;
U_(i,t,iter)=U.l(i,t);
S_(pv,t,iter)=S.l(pv,t);
lowerbound=z_master.l;
);


if(s_optcut.l=1 and s_fcut.l=0,
solve master_problem_optcut minimizing z_master using MIP;
display z_master.l,U.l, S.l, Pg.l;
U_(i,t,iter)=U.l(i,t);
S_(pv,t,iter)=S.l(pv,t);
lowerbound=z_master.l;
);
if (s_fcut.l=1 and s_optcut.l=0,
solve master_problem_fcut minimizing z_master using MIP;
display z_master.l,U.l,S.l, Pg.l;
U_(i,t,iter)=U.l(i,t);
S_(pv,t,iter)=S.l(pv,t);
lowerbound=z_master.l;
);
if (s_optcut.l=1 and s_fcut.l=1,
solve master_problem minimizing z_master using MIP;
display z_master.l,U.l, S.l;
U_(i,t,iter)=U.l(i,t);
S_(pv,t,iter)=S.l(pv,t);
lowerbound=z_master.l;
);


option LP=CPLEX;
option optcr=0;
option threads=4;
solve feasiblity_check minimizing z_ckeck using LP;
vt_(t,iter)=vt.l(t);

if(LS_sum.l=0,
option LP=CPLEX;
option threads=4;
solve  optimal_operation minimizing z_opt using lp;
optcut(iter)=yes;
s_optcut.l=1;
uperbound=z_opt.l;
pma_(pv,t,iter)=pma.m(pv,t);
pmi_(pv,t,iter)=pmi.m(pv,t);

else
fcut(iter)=yes;
s_fcut.l=1;
pma_check_(pv,t,iter)=pma_check.m(pv,t);
pmi_check_(pv,t,iter)=pmi_check.m(pv,t);
);
results(iter,'lb')=lowerbound;
results(iter,'ub')=uperbound;
nonconverged$((uperbound-lowerbound)<=0)=no;
);

display  uperbound, lowerbound, Power_Flow.m;

execute_unload 'UC_SFR problem.gdx'


*==============================================================================================================================================================================================================================================================================
*execute_unload "UC_SFRn30.gdx" z.l
*execute  'gdxxrw.exe UC_SFRn30.gdx  var=z  rng=z!a1'
*execute_unload "UC_SFRn30.gdx" U.l
*execute  'gdxxrw.exe UC_SFRn30.gdx  var=U  rng=U!a1'
execute_unload "UC_SFRn30.gdx" Power_Flow.m
execute  'gdxxrw.exe UC_SFRn30.gdx  var=Power_Flow.m  rng=Power_Flow.m!a1'
*execute_unload "UC_SFRn30.gdx" f.l
*execute  'gdxxrw.exe UC_SFRn30.gdx  var=f  rng=f!a1'

*==============================================================================================================================================================================================================================================================================
*Display zdsp.l,pg_pv_te.m,Pg_pv_t_end.m,teta.m;
