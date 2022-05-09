**** OptKnock: 	A Bilevel Programming Framework for Identifying Gene Knockout Strategies for 
****		Microbial Strain Optimization

**** Original Author(s): 	Anthony P. Burgard, Priti Pharkya and Costas D. Maranas
**** Developed at the 		Chemical & Biological Systems Optimizations Laboratory
****				The Pennsylvania State University, University Park, PA, USA
**** Code Redeveloped by: 	Sridhar Ranganathan, The Pennsylvania State University
**** Contact: 			costas@engr.psu.edu, sur152@psu.edu



**** Notes:
*    	1. Run this file after estimating the initial bounds
*    	2. The initial bounds are set as parameters - basemin and basemax		

$set myroot iAF1260/iAF1260

options
        limrow = 10000
        limcol = 10000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
        sysout = off
        solprint = on
        mip = OSIGUROBI;

sets 

iter /1*25/

i metabolites 
$include "%myroot%_cmp.txt"

j reactions 
$include "%myroot%_rxnnames.txt"

$ONEMPTY
sourcemetab(i)
$include "%myroot%_source.txt"

sourceflux(j)
$include "%myroot%_sourceflux.txt"

offaero(j)
$include "%myroot%_offaero.txt"

offglu(j)
$include "%myroot%_offglu.txt"

offaeroglu(j)
$include "%myroot%_offaeroglu.txt"

fwdrev(j) 

constraints(j)
/
'EX_glc(e)'
'Ec_biomass_iAF1260_WT_59p81M'
/

spontaneous(j)
/
'ACALDtpp_f'
'ACONIs_f'
'AOBUTDs'
'ARBTNexs'
'ATPHs'
'CO2tpp_f'
'CPGNexs'
'DATPHs'
'DHPTDCs'
'FALDtpp_f'
'FALGTHLs_f'
'FE3HOXexs'
'FECRMexs'
'FEENTERexs'
'FEOXAMexs'
'G5SADs'
'GLYCtpp_f'
'GTPHs'
'H2Otex_f'
'H2Otpp_f'
'H2St1pp'
'H2tpp_f'
'METOX1s'
'METOX2s'
'N2Otpp_f'
'NH4tpp_f'
'NOtpp_f'
'O2tpp_f'
'SO2tpp_f'
'SUCCtex_f'
'EX_hdca(e)'
/

transporters(j)
$include "%myroot%_transporters.txt"

;

alias(j, j1);
alias(iter, iter1);

parameter 

S(i,j)
$include "%myroot%_sij.txt"

Sneg(i, j)

rxntype(j) reaction type
$include "%myroot%_rxntype.txt" 

basemin(j)
$include "%myroot%Min.txt"

basemax(j)
$include "%myroot%Max.txt"

LB(j), UB(j)

k
temp(j)
blocked(j)
countstop
objdiff
bioprod
knockoutsyet
count
ofvalue(iter)
allowknock
notstop
objstore(iter)
ystore(iter, j)
wildtype
*y(j)
;

scalar minbiomass /5/;
scalar M /1000/;
*y(j) = 1;
k = 0;
Sneg(i, j) = -S(i, j);

fwdrev(j)$(rxntype(j) = 1 or rxntype(j) = 2) = yes;

variables 
v(j)
zprimal
zdual
z
mu(j)
lambda(i)
alpha1, beta1
;

positive variables
alpha(j)
beta(j)
biomass
;

binary variables y(j);

equations

primalobj
primal1
primal2
primal3
primal4
primal5

dualobj
dual1
dual2
dual3

outerobj
outer1
outer2
outer3
outer4
outer5
outer6
outer7
outer8
outer9
outer10
outer11
;

*** PRIMAL INNER PROBLEM 
*** MAXIMIZE BIOMASSS

primalobj..					zprimal =e= v('Ec_biomass_iAF1260_WT_59p81M');
primal1(i)$(not sourcemetab(i))..		sum(j, (S(i,j)*v(j))) =e= 0;
primal2..					v('EX_glc(e)') =e= -100;
primal3..					v('Ec_biomass_iAF1260_WT_59p81M') =g= minbiomass;
primal4(j)..					v(j) =l= basemax(j)*y(j);
primal5(j)..					-v(j) =l= -basemin(j)*y(j);

*** DUAL INNER PROBLEM

dualobj..		zdual =e= - mu('EX_glc(e)')*100 - minbiomass*biomass + sum(j, alpha(j)*basemax(j)*y(j) - beta(j)*basemin(j)*y(j) );
dual1..			sum(i$(not sourcemetab(i)), lambda(i)*S(i,'Ec_biomass_iAF1260_WT_59p81M') ) - biomass + alpha('Ec_biomass_iAF1260_WT_59p81M') - beta('Ec_biomass_iAF1260_WT_59p81M') =e= 1;
dual2..			sum(i$(not sourcemetab(i)), lambda(i)*S(i,'EX_glc(e)') ) + mu('EX_glc(e)') + alpha('EX_glc(e)') - beta('EX_glc(e)') =e= 0;
dual3(j)$(not constraints(j))..		
			sum(i$(not sourcemetab(i)), lambda(i)*S(i, j) ) + alpha(j) - beta(j) =e= 0;

*** OUTER PROBLEM
*** MAXIMIZE V(succinate)

outerobj..		z =e= v('THRS');
outer1..		sum(j, (1-y(j) ) ) =l= allowknock;
outer2..		v('Ec_biomass_iAF1260_WT_59p81M') =e= - mu('EX_glc(e)')*100 - minbiomass*biomass + sum(j, alpha1(j)*basemax(j) - beta1(j)*basemin(j) );

*** Linearinzing constraints

outer3(j)..		alpha1(j) =l= M*y(j);
outer4(j)..		alpha1(j) =g= -M*y(j);
outer5(j)..		alpha1(j) =l= alpha(j) + M*(1-y(j));
outer6(j)..		alpha1(j) =g= alpha(j) - M*(1-y(j));

outer7(j)..             beta1(j) =l= M*y(j);
outer8(j)..             beta1(j) =g= -M*y(j);
outer9(j)..             beta1(j) =l= beta(j) + M*(1-y(j));
outer10(j)..            beta1(j) =g= beta(j) - M*(1-y(j));

*** Integer Cut Constraint

outer11(iter1)$((ord(iter1) lt count) and knockoutsyet)..	v('THRS') =g= objstore(iter1) + 0.01 - 1000*sum(j$(ystore(iter1, j) lt 0.5), y(j));

*** Lower and Upper Bounds for the fluxes

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = 1000;
LB(j)$(rxntype(j) = 1) = -1000;
UB(j)$(rxntype(j) = 1) = 1000;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = 1000;
LB(j)$(rxntype(j) = 3) = -1000;
UB(j)$(rxntype(j) = 3) = 1000;

LB('EX_glc(e)') = -100;
LB('EX_o2(e)') = -1000;

*anaerobic glucose uptake conditions
LB(j)$(offglu(j)) = 0;
UB(j)$(offglu(j)) = 0;

v.lo(j) = LB(j);
v.up(j) = UB(j);
y.fx(j)$(rxntype(j) = 2) = 1;
y.fx(j)$(spontaneous(j)) = 1;
y.fx(j)$(transporters(j)) = 1;

model primalproblem
/
primal1
primal2
primal3
primal4
primal5
primalobj
/
;

model dualproblem 
/
dual1
dual2
dual3
dualobj
/
;

model optknock
/
primal1
primal2
primal3
primal4
primal5
dual1
dual2
dual3
outerobj
outer1
outer2
outer3
outer4
outer5
outer6
outer7
outer8
outer9
outer10
outer11
/
;

*solve primalproblem using lp maximizing zprimal;

*solve dualproblem using lp minimizing zdual;

file res /OptKnock_Results.txt/;

blocked(j)$( (UB(j) eq 0) and (LB(j) eq 0) ) = 1;
S(i, j)$(blocked(j)) = 0;

y.fx(j)$( (LB(j) gt 0) and (UB(j) gt 0) ) = 1;
y.fx(j)$( (LB(j) lt 0) and (UB(j) lt 0) ) = 1;
y.fx(j)$(blocked(j)) = 1;
y.fx(j)$(rxntype(j) = 4) = 1;
y.fx(j)$(rxntype(j) = 2) = 1;
y.fx(j)$(rxntype(j) = 3) = 1;

ystore(iter, j) = 0;
count = 0;
notstop = 1;
allowknock = 3;
ofvalue(iter) = 0;
objstore(iter) = 0;
knockoutsyet = 0;
wildtype = 0;

loop(iter$notstop,
       count = ord(iter);

       y.l(j) = 1;
       v.l(j) = 0;

       optknock.optfile = 1;

       solve optknock using mip maximizing z ;

       ofvalue(iter) = 1;

       if ( (knockoutsyet eq 0),
         wildtype = z.l;
       );

       objstore(iter) = v.l('THRS');
       knockoutsyet = 1;
       objdiff = v.l('THRS') - wildtype;


       if((optknock.modelstat ne 1),
               ystore(iter,j) = 1;
               countstop = count + 1;  
               notstop = 0;);

       if( (objdiff lt 10e-3),
               ofvalue(iter) = 0;
               allowknock = allowknock + 1);

       if ( (allowknock gt 7),
               countstop = count + 1;
               notstop = 0;);

       if(((optknock.modelstat eq 1) and ((objdiff gt 10e-3) or (knockoutsyet eq 0)) ),
               ofvalue(iter) = 1;
               ystore(iter,j)$(y.l(j) gt 0.5) = 1;
               ystore(iter,j)$(y.l(j) lt 0.5) = 0;);

       put res;

 PUT "Threonine  = ",z.l:12:6/;
 put "Biomass    = ",v.l('Ec_biomass_iAF1260_WT_59p81M'):12:6/;
 Put "Model Status =", optknock.modelstat:2:0 /;
 put "CPU for MIP:", optknock.resusd:6:4 ///;

* LOOP(j$(v.l(j) ne 0),
* *PUT "v(",j.tl:7:0,") = ",v.l(j):10:5," ";
* *  if( (v.l(j) gt 0),
* *      loop(i$(S(i,j) lt 0),
*         if ( (Sneg(i,j) eq 1),        
*           put i.tl:6:0," ";)
*         if ( (Sneg(i,j) ne 1),
*           put Sneg(i,j):5:2," ",i.tl:7:0," ";););    
*       put "-->";
*       loop(i$(S(i,j) gt 0),
*         if ( (S(i,j) eq 1),
*           put i.tl:6:0," ";)
*         if ( (S(i,j) ne 1),
*         put S(i,j):5:2," ",i.tl:7:0," ";););
*   );
*    if( (v.l(j) lt 0),
*       loop(i$(S(i,j) gt 0),
*         if ( (S(i,j) eq 1),
*           put i.tl:6:0," ";)
*         if ( (S(i,j) ne 1),
*           put S(i,j):5:2," ",i.tl:7:0," ";););
*       put "-->";
*       loop(i$(S(i,j) lt 0),
*         if ( (Sneg(i,j) eq 1),
*           put i.tl:6:0," ";)
*         if ( (Sneg(i,j) ne 1),
*           put Sneg(i,j):5:2," ",i.tl:7:0," ";););
*   );
*PUT /;
* );

loop(j$(y.l(j) eq 0),
       put "y(",j.tl:10:0,")= zero"/;);
put "******************************************"/;
put /;


);

loop(iter$ofvalue(iter),
 put "Objstore(",iter.tl:4:0,")= ",objstore(iter):10:5/;
 loop(j$(ystore(iter,j) eq 0),
       put "  y(",j.tl:10:0,")= ",ystore(iter,j):4:0,"     ";
       loop(i$(S(i,j) lt 0),
         if ( (Sneg(i,j) eq 1),
           put i.tl:6:0," ";)
         if ( (Sneg(i,j) ne 1),
           put Sneg(i,j):5:2," ",i.tl:7:0," ";););
       put "-->";
       loop(i$(S(i,j) gt 0),
         if ( (S(i,j) eq 1),
           put i.tl:6:0," ";)
         if ( (S(i,j) ne 1),
         put S(i,j):5:2," ",i.tl:7:0," ";););
       put /;
 );
 put /;);
