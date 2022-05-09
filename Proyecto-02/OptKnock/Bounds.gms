**** OptKnock code **************************************
* 	Notes:
*	This code finds the initial bounds for the 
*	metabolic network given the production of
*	biomass is greater than a prespecified value.

$set myroot iAF1260/iAF1260

options limrow = 1000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
        mip = OSIGUROBI;

sets 
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
vj(j)
;

parameter 

* stoichiometric coefficients
s(i,j)
$include "%myroot%_sij.txt"

* reaction type
rxntype(j)
$include "%myroot%_rxntype.txt"

n  
basemin(j)
basemax(j)
;

fwdrev(j)$(rxntype(j) = 1 or rxntype(j) = 2) = yes;

variables v(j), z;

equations 

network, obj
primal1, primal2
;

network(i)$(not sourcemetab(i))..	sum(j, (s(i,j)*v(j))) =e= 0;
obj..                                   z =e= (sum(j$(vj(j)), v(j)));
primal1..                               v('EX_glc(e)') =e= -100;
primal2..                               v('Ec_biomass_iAF1260_WT_59p81M ') =g= 5;

parameter LB(j), UB(j);
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

model flux_base
/
all
/
;

flux_base.optfile=1;

for (n=1 to card(j) by 1,

v.lo(j)=LB(j);
v.up(j)=UB(j);

vj(j)$(ord(j) eq n) = yes;
vj(j)$(ord(j) ne n) = no;
solve flux_base using lp minimizing z;
basemin(vj) = z.l;
solve flux_base using lp maximizing z;
basemax(vj) = z.l;
);

file base1 /Bounds.txt/;
put base1;

loop(j,
if (rxntype(j) ne 2,
put j.tl, basemin(j):15:8, basemax(j):15:8,;
if (basemin(j) eq basemax(j),
put '   FIXED'/;
elseif (basemin(j) eq -1000 and basemax(j) eq 1000),
put '   LIMITS'/;
else
put /;
);
);
);

file base_min /iAF1260Min.txt/;
put base_min;
put "/"/;
loop(j,
put "'"j.tl:30:0"'", basemin(j):9:5/;
);
put "/";

file base_max /iAF1260Max.txt/;
put base_max;
put '/'/;
loop(j,
put "'"j.tl:30:0"'", basemax(j):9:5/;
);
put "/";
