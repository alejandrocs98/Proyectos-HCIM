$set myroot iAF1260/iAF1260

options limrow = 1000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
        mip = OSIGUROBI;

sets i metabolites
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

fwdrev(j), 
vj(j)
;

parameter s(i,j)
$include "%myroot%_sij.txt"

rxntype(j) reaction type
$include "%myroot%_rxntype.txt"

n 
basemin(j), basemax(j)
vmin(j), vmax(j)

fixed, limit, range
;

fwdrev(j)$(rxntype(j) = 1 or rxntype(j) = 2) = yes;

variables v(j), z;

equations 

network, obj
glc_upt 
; 

network(i)$(not sourcemetab(i))..	sum(j, (s(i,j)*v(j))) =e= 0;
obj..					z =e= v('Ec_biomass_iAF1260_core_59p81M');
glc_upt..				v('EX_glc(e)') =e= -100;
	
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

solve flux_base using lp maximizing z;
