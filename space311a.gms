SET	t	 time periods	/0*100/;

SETS
	tfirst(t)     first period
	tlast(t)      last period
    z             location  /earth, space/
    scn           /1*7/
    szn           /9*10/
;

tfirst(t) = yes$(ord(t) eq 1);
tlast(t) = yes$(ord(t) eq card(t));

PARAMETERS
	g	population growth rate				        /0.023/
	delta	depreciation rate 				        /0.04/
    deltav  knowledge depreciation                  /0.0095/
	i0(z)	initial investment				        /earth 0.3, space 1e-4/
	c0	initial consumption				            /0.22/
	eta	inverse of the intertemporal elasticity		/2/

*Production Additions:
	r0(z)  initial RD						        /earth 0.1, space 1e-4/
*	b(z)	value share of capital in C-D           /earth 0.5, space 0.7/
*SA #1:
	b(z)	value share of capital in C-D           /earth 0.5, space 0.7/
*S-Curve Parameters:
	sig(z)	upper asymptote of scalar		        /earth 0.7, space 1.3/
	nu	divide upper asymptote			            /1/
	chi(z)	horizontal stretch of scalar	        /earth 1, space 0.95/
	zeta	lower asymptote of scalar
	alpha(z)	middel of S-curve
*CO2 Parameters:
*    carblim max carbon stock                        /10/
    carblim max carbon stock                        
    e(z)    emissions intensity                     /earth 0.1, space 0/
    ds  CO2 dissipation rate                        /0.05/
;

PARAMETERS
	rho	    pure time preference rate
	k0(z)	initial capital stock
	V0(z)	initial knowledge stock
	l0	    initial labor supply
	a0(z)      initial scalar parameter in C-D prod fn
	qref(t) steady-state output index
	
	l(t)	labor supply
	beta(t)	social welfare weight

    ci(z)   capital investment cost /earth 0.55, space 1.9/
    cr(z)   r&D cost                /earth 2, space 2/

	a1		space production share	/0.3/
	a2		earth production share	/0.7/
	psi		consumption good elasticity of sub	/0.5/

*Sensitivity analyses
    init    initial value for solves
    SA1     b_space     /0.8/
    SA2     b_space     /0.6/
    SA3     b_earth     /0.6/
    SA4     b_earth     /0.4/
    SA5     ci_space    /2.5/
    SA6     ci_space    /1.1/
    SA7     ci_earth    /0.9/
    SA8     ci_earth    /0.3/

    capset(scn)  alt carbon ceilings /1 5, 2 6, 3 7, 4 8, 5 9, 6 10, 7 11/
    alphalt(szn) alt s-curve space spec /9 5.7, 10 5.9/
;

$if not setglobal var1 $setglobal var1 0
$if not setglobal var1 $setglobal var2 0

carblim = %var1%;

display carblim;

alpha("space") = %var2%;

display alpha;

alpha("earth") = 4.0;

k0(z) = i0(z)/(g+delta);
V0(z) = r0(z)/(g+deltav);

qref(t) = power(1+g, ord(t)-1);

l0 = (1-b("earth"))*(c0+sum(z,i0(z)+r0(z)));
l(t) = l0 * qref(t);

zeta = (c0 + sum(z, i0(z) + r0(z)))/sum(z, k0(z)**b("earth")*l0**(1-b("earth"))) - sig("earth")/(nu+ exp(-chi("earth") * sum(z, V0(z)) + alpha("earth")));

a0(z) = zeta + sig(z) / (nu + exp(-chi(z) * V0(z) + alpha(z)));

*xx change this to so b is summed over z
rho = a0("earth") * b("earth") * sum(z,k0(z))**(b("earth")-1) * l0**(1-b("earth")) - delta;

abort$(g > rho) "Growth rate exceeds the discount rate!!!", g, rho;

beta(t) = power((1+g)**eta/(1+rho), ord(t)-1);

VARIABLES
	C(t)	consumption
	K(z,t)	capital stock
	Y(z,t)	output
	I(z,t)	investment
	W	social welfare

	V(z,t)	knowledge stocks
	R(z,t)	RD investment
    a(z,t)    production scalar
    m(t)    atmospheric carbon level
    lz(z,t) location-specific labor

	pearth(t)	earth output used for consumption good
	pspace(t)	space output used for consumption good
    it(t)       total investment
    rt(t)       total R&D
;

EQUATIONS
	production(z,t)	definition of prod fn
	allocation(t)	output accounting identity
	accumulation(z,t)	capital stock evolution
	utility		social welfare fn

	vstock(z,t)		knowledge eqn of motion
    scal(z,t)         production scalar eqn

    atm(t)          atmospheric carbon level eqn
    ceil(t)         carbon ceiling eqn
    lab(t)          labor allocation

	cons(t)			consumption good production fn
    itot(t)         total investment
    rtot(t)         total R&D
;

utility..		    W =e= sum(t, beta(t) * C(t)**(1-eta) / (1-eta) );

production(z,t).. 	Y(z,t) =e= a(z,t) * K(z,t)**b(z) * lz(z,t)**(1-b(z));

cons(t)..			C(t) =e= (a1* pspace(t)**psi + a2* pearth(t)*psi)**(1/psi);

allocation(t)..		pspace(t) + pearth(t) + sum(z, ci(z)* I(z,t) + cr(z)* R(z,t)) =l= sum(z, Y(z,t));

accumulation(z,t)..	K(z,t) =e= (1-delta)*K(z,t-1) + I(z,t-1) + k0(z)$tfirst(t);

vstock(z,t)..		V(z,t) =e= (1-deltav)*V(z,t-1) + R(z,t-1) + V0(z)$tfirst(t);

scal(z,t)..         a(z,t) =e= zeta + sig(z) / (nu + exp(-chi(z) * V(z,t) + alpha(z)));

atm(t)..            m(t) =e= m(t-1) - ds + sum(z, Y(z,t)*e(z));

ceil(t)..           m(t) =l= carblim;

lab(t)..            sum(z, lz(z,t)) =l= l(t);

itot(t)..           it(t) =e= sum(z, I(z,t));
rtot(t)..           rt(t) =e= sum(z, R(z,t));

C.l(t) = c0 * qref(t);
I.l("earth",t) = i0("earth") * qref(t);
it.l(t) = i0("earth")*qref(t);
rt.l(t) = r0("earth")*qref(t);
Y.l("earth",t) = (c0+i0("earth")+r0("earth")) * qref(t);

R.l("earth",t) = R0("earth") * qref(t);

C.lo(t) = 1e-4;
I.lo(z,t) = 1e-4;
K.lo(z,t) = 1e-4;

R.lo(z,t) = 1e-4;
V.lo(z,t) = 1e-4;
lz.lo(z,t) = 1e-4;
Y.lo(z,t) = 0;

pearth.lo(t) = 1e-4;
pspace.lo(t) = 1e-4;

MODEL	 noterm	   /utility, production, allocation, accumulation, vstock, scal, atm, ceil, lab, cons, itot, rtot/;

solve noterm maximizing W using NLP;

    display alpha;
    display a0;
    display a.l;

PARAMETERS
	report output with t and x index
	reportt output with t index
;

report(t, z, "K", "%var1%", "%var2%") = K.l(z,t);
report(t, z, "R", "%var1%", "%var2%") = R.l(z,t);
report(t, z, "V", "%var1%", "%var2%") = V.l(z,t);
report(t, z, "Y", "%var1%", "%var2%") = Y.l(z,t);
report(t, z, "I", "%var1%", "%var2%") = I.l(z,t);
report(t, z, "a", "%var1%", "%var2%") = a.l(z,t);
report(t, z, "C", "%var1%", "%var2%") = C.l(t);
report(t, z, "alpha", "%var1%", "%var2%") = alpha(z);
report(t, z, "m", "%var1%", "%var2%") = m.l(t);
report(t, z, "ko(z)", "%var1%", "%var2%") = k0(z);
report(t, z, "v0(z)", "%var1%", "%var2%") = v0(z);
report(t, z, "lz(z,t)", "%var1%", "%var2%") = lz.l(z,t);

reportt(t, "rho", "%var1%", "%var2%") = rho;
reportt(t, "qref", "%var1%", "%var2%") = qref(t);
reportt(t, "l0", "%var1%", "%var2%") = rho;
reportt(t, "lt", "%var1%", "%var2%") = l(t);
reportt(t, "zeta", "%var1%", "%var2%") = zeta;

execute_unload 'baseramsey_%var1%_%var2%.gdx' report, reportt;
$exit
alpha("space") = 5.7;
a0(z) = zeta + sig(z) / (nu + exp(-chi(z) * V0(z) + alpha(z)));
MODEL	 TESTM	   /utility, production, allocation, accumulation, vstock, scal, atm, ceil, lab, cons, itot, rtot/;

solve TESTM maximizing W using NLP;

    display alpha;
    display a0;
    display a.l;

PARAMETERS
	report output with t and x index
	reportt output with t index
;

report(t, z, "K", "test") = K.l(z,t);
report(t, z, "R", "test") = R.l(z,t);
report(t, z, "V", "test") = V.l(z,t);
report(t, z, "Y", "test") = Y.l(z,t);
report(t, z, "I", "test") = I.l(z,t);
report(t, z, "a", "test") = a.l(z,t);
report(t, z, "C", "test") = C.l(t);
report(t, z, "alpha", "test") = alpha(z);
report(t, z, "m", "test") = m.l(t);
report(t, z, "ko(z)", "test") = k0(z);
report(t, z, "v0(z)", "test") = v0(z);
report(t, z, "lz(z,t)", "test") = lz.l(z,t);

reportt(t, "rho", "test") = rho;
reportt(t, "qref", "test") = qref(t);
reportt(t, "l0", "test") = rho;
reportt(t, "lt", "test") = l(t);
reportt(t, "zeta", "test") = zeta;

init = alpha("space");

loop(szn,
    alpha("space") = alphalt(szn);
    solve noterm maximizing W using NLP;
    display alpha;
    display a0;
    display a.l;

    report(t, z, "K", szn) = K.l(z,t);
    report(t, z, "R", szn) = R.l(z,t);
    report(t, z, "V", szn) = V.l(z,t);
    report(t, z, "Y", szn) = Y.l(z,t);
    report(t, z, "I", szn) = I.l(z,t);
    report(t, z, "a", szn) = a.l(z,t);
    report(t, z, "C", szn) = C.l(t);
    report(t, z, "alpha", szn) = alpha(z);
    report(t, z, "m", szn) = m.l(t);
    report(t, z, "ko(z)", szn) = k0(z);
    report(t, z, "v0(z)", szn) = v0(z);
    report(t, z, "lz(z,t)", szn) = lz.l(z,t);
    
    reportt(t, "rho", szn) = rho;
    reportt(t, "qref", szn) = qref(t);
    reportt(t, "l0", szn) = rho;
    reportt(t, "lt", szn) = l(t);
    reportt(t, "zeta", szn) = zeta;
);
*$offtext

alpha("space") = init;

execute_unload "baseramsey.gdx" report, reportt;
$exit

*SA1:
init = b("space");
b("space") = SA1;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA1") = K.l(z,t);
report(t, z, "R", "SA1") = R.l(z,t);
report(t, z, "V", "SA1") = V.l(z,t);
report(t, z, "Y", "SA1") = Y.l(z,t);
report(t, z, "I", "SA1") = I.l(z,t);
report(t, z, "a", "SA1") = a.l(z,t);
report(t, z, "C", "SA1") = C.l(t);
report(t, z, "alpha", "SA1") = alpha(z);
report(t, z, "m", "SA1") = m.l(t);
report(t, z, "ko(z)", "SA1") = k0(z);
report(t, z, "v0(z)", "SA1") = v0(z);
report(t, z, "lz(z,t)", "SA1") = lz.l(z,t);

reportt(t, "rho", "SA1") = rho;
reportt(t, "qref", "SA1") = qref(t);
reportt(t, "l0", "SA1") = rho;
reportt(t, "lt", "SA1") = l(t);
reportt(t, "zeta", "SA1") = zeta;

b("space") = init;

*SA2:
init = b("space");
b("space") = SA2;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA2") = K.l(z,t);
report(t, z, "R", "SA2") = R.l(z,t);
report(t, z, "V", "SA2") = V.l(z,t);
report(t, z, "Y", "SA2") = Y.l(z,t);
report(t, z, "I", "SA2") = I.l(z,t);
report(t, z, "a", "SA2") = a.l(z,t);
report(t, z, "C", "SA2") = C.l(t);
report(t, z, "alpha", "SA2") = alpha(z);
report(t, z, "m", "SA2") = m.l(t);
report(t, z, "ko(z)", "SA2") = k0(z);
report(t, z, "v0(z)", "SA2") = v0(z);
report(t, z, "lz(z,t)", "SA2") = lz.l(z,t);

reportt(t, "rho", "SA2") = rho;
reportt(t, "qref", "SA2") = qref(t);
reportt(t, "l0", "SA2") = rho;
reportt(t, "lt", "SA2") = l(t);
reportt(t, "zeta", "SA2") = zeta;

b("space") = init;

*SA3:
init = b("earth");
b("earth") = SA3;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA3") = K.l(z,t);
report(t, z, "R", "SA3") = R.l(z,t);
report(t, z, "V", "SA3") = V.l(z,t);
report(t, z, "Y", "SA3") = Y.l(z,t);
report(t, z, "I", "SA3") = I.l(z,t);
report(t, z, "a", "SA3") = a.l(z,t);
report(t, z, "C", "SA3") = C.l(t);
report(t, z, "alpha", "SA3") = alpha(z);
report(t, z, "m", "SA3") = m.l(t);
report(t, z, "ko(z)", "SA3") = k0(z);
report(t, z, "v0(z)", "SA3") = v0(z);
report(t, z, "lz(z,t)", "SA3") = lz.l(z,t);

reportt(t, "rho", "SA3") = rho;
reportt(t, "qref", "SA3") = qref(t);
reportt(t, "l0", "SA3") = rho;
reportt(t, "lt", "SA3") = l(t);
reportt(t, "zeta", "SA3") = zeta;

b("earth") = init;

*SA4:
init = b("earth");
b("earth") = SA4;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA4") = K.l(z,t);
report(t, z, "R", "SA4") = R.l(z,t);
report(t, z, "V", "SA4") = V.l(z,t);
report(t, z, "Y", "SA4") = Y.l(z,t);
report(t, z, "I", "SA4") = I.l(z,t);
report(t, z, "a", "SA4") = a.l(z,t);
report(t, z, "C", "SA4") = C.l(t);
report(t, z, "alpha", "SA4") = alpha(z);
report(t, z, "m", "SA4") = m.l(t);
report(t, z, "ko(z)", "SA4") = k0(z);
report(t, z, "v0(z)", "SA4") = v0(z);
report(t, z, "lz(z,t)", "SA4") = lz.l(z,t);

reportt(t, "rho", "SA4") = rho;
reportt(t, "qref", "SA4") = qref(t);
reportt(t, "l0", "SA4") = rho;
reportt(t, "lt", "SA4") = l(t);
reportt(t, "zeta", "SA4") = zeta;

b("earth") = init;

*SA5:
init = ci("space");
ci("space") = SA5;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA5") = K.l(z,t);
report(t, z, "R", "SA5") = R.l(z,t);
report(t, z, "V", "SA5") = V.l(z,t);
report(t, z, "Y", "SA5") = Y.l(z,t);
report(t, z, "I", "SA5") = I.l(z,t);
report(t, z, "a", "SA5") = a.l(z,t);
report(t, z, "C", "SA5") = C.l(t);
report(t, z, "alpha", "SA5") = alpha(z);
report(t, z, "m", "SA5") = m.l(t);
report(t, z, "ko(z)", "SA5") = k0(z);
report(t, z, "v0(z)", "SA5") = v0(z);
report(t, z, "lz(z,t)", "SA5") = lz.l(z,t);

reportt(t, "rho", "SA5") = rho;
reportt(t, "qref", "SA5") = qref(t);
reportt(t, "l0", "SA5") = rho;
reportt(t, "lt", "SA5") = l(t);
reportt(t, "zeta", "SA5") = zeta;

ci("space") = init;

*SA6:
init = ci("space");
ci("space") = SA6;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA6") = K.l(z,t);
report(t, z, "R", "SA6") = R.l(z,t);
report(t, z, "V", "SA6") = V.l(z,t);
report(t, z, "Y", "SA6") = Y.l(z,t);
report(t, z, "I", "SA6") = I.l(z,t);
report(t, z, "a", "SA6") = a.l(z,t);
report(t, z, "C", "SA6") = C.l(t);
report(t, z, "alpha", "SA6") = alpha(z);
report(t, z, "m", "SA6") = m.l(t);
report(t, z, "ko(z)", "SA6") = k0(z);
report(t, z, "v0(z)", "SA6") = v0(z);
report(t, z, "lz(z,t)", "SA6") = lz.l(z,t);

reportt(t, "rho", "SA6") = rho;
reportt(t, "qref", "SA6") = qref(t);
reportt(t, "l0", "SA6") = rho;
reportt(t, "lt", "SA6") = l(t);
reportt(t, "zeta", "SA6") = zeta;

ci("space") = init;

*SA7:
init = ci("earth");
ci("earth") = SA7;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA7") = K.l(z,t);
report(t, z, "R", "SA7") = R.l(z,t);
report(t, z, "V", "SA7") = V.l(z,t);
report(t, z, "Y", "SA7") = Y.l(z,t);
report(t, z, "I", "SA7") = I.l(z,t);
report(t, z, "a", "SA7") = a.l(z,t);
report(t, z, "C", "SA7") = C.l(t);
report(t, z, "alpha", "SA7") = alpha(z);
report(t, z, "m", "SA7") = m.l(t);
report(t, z, "ko(z)", "SA7") = k0(z);
report(t, z, "v0(z)", "SA7") = v0(z);
report(t, z, "lz(z,t)", "SA7") = lz.l(z,t);

reportt(t, "rho", "SA7") = rho;
reportt(t, "qref", "SA7") = qref(t);
reportt(t, "l0", "SA7") = rho;
reportt(t, "lt", "SA7") = l(t);
reportt(t, "zeta", "SA7") = zeta;

ci("earth") = init;

*SA8:
init = ci("earth");
ci("earth") = SA8;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA8") = K.l(z,t);
report(t, z, "R", "SA8") = R.l(z,t);
report(t, z, "V", "SA8") = V.l(z,t);
report(t, z, "Y", "SA8") = Y.l(z,t);
report(t, z, "I", "SA8") = I.l(z,t);
report(t, z, "a", "SA8") = a.l(z,t);
report(t, z, "C", "SA8") = C.l(t);
report(t, z, "alpha", "SA8") = alpha(z);
report(t, z, "m", "SA8") = m.l(t);
report(t, z, "ko(z)", "SA8") = k0(z);
report(t, z, "v0(z)", "SA8") = v0(z);
report(t, z, "lz(z,t)", "SA8") = lz.l(z,t);

reportt(t, "rho", "SA8") = rho;
reportt(t, "qref", "SA8") = qref(t);
reportt(t, "l0", "SA8") = rho;
reportt(t, "lt", "SA8") = l(t);
reportt(t, "zeta", "SA8") = zeta;

ci("earth") = init;

execute_unload "baseramsey.gdx" report, reportt;