SET	t	 time periods	/0*100/;

SETS
	tfirst(t)     first period
	tlast(t)      last period
    z             location  /earth, space/
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
*SA #1:
	b(z)	value share of capital in C-D           /earth 0.5, space 0.7/
*S-Curve Parameters:
	sig(z)	upper asymptote of scalar		        /earth 0.7, space 1.3/
	nu	divide upper asymptote			            /1/
	chi(z)	horizontal stretch of scalar	        /earth 1, space 0.95/
	zeta	lower asymptote of scalar
	alpha(z)	middel of S-curve                   /earth 4, space 5.1/
*CO2 Parameters:
    carblim max carbon stock                        /10/
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

	a1		earth production share	/0.7/
	a2		space production share	/0.3/
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
	SA9     cr_space    /3/
    SA10    cr_space    /1/
    SA11    cr_earth    /3/
    SA12    cr_earth    /1/
;

k0(z) = i0(z)/(g+delta);
V0(z) = r0(z)/(g+deltav);

qref(t) = power(1+g, ord(t)-1);

l0 = (1-b("earth"))*(c0+sum(z,i0(z)+r0(z)));
l(t) = l0 * qref(t);

zeta = (c0 + sum(z, i0(z) + r0(z)))/sum(z, k0(z)**b("earth")*l0**(1-b("earth"))) - sig("earth")/(nu+ exp(-chi("earth") * sum(z, V0(z)) + alpha("earth")));

a0(z) = zeta + sig(z) / (nu + exp(-chi(z) * V0(z) + alpha(z)));

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

PARAMETERS
	report output with t and x index
	reportt output with t index
;

report(t, z, "K", "Base") = K.l(z,t);
report(t, z, "R", "Base") = R.l(z,t);
report(t, z, "V", "Base") = V.l(z,t);
report(t, z, "Y", "Base") = Y.l(z,t);
report(t, z, "I", "Base") = I.l(z,t);
report(t, z, "a", "Base") = a.l(z,t);
report(t, z, "C", "Base") = C.l(t);
report(t, z, "alpha", "Base") = alpha(z);
report(t, z, "m", "Base") = m.l(t);
report(t, z, "ko(z)", "Base") = k0(z);
report(t, z, "v0(z)", "Base") = v0(z);
report(t, z, "lz(z,t)", "Base") = lz.l(z,t);

reportt(t, "rho", "Base") = rho;
reportt(t, "qref", "Base") = qref(t);
reportt(t, "l0", "Base") = rho;
reportt(t, "lt", "Base") = l(t);
reportt(t, "zeta", "Base") = zeta;

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

*SA9
init = cr("space");
cr("space") = SA9;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA9") = K.l(z,t);
report(t, z, "R", "SA9") = R.l(z,t);
report(t, z, "V", "SA9") = V.l(z,t);
report(t, z, "Y", "SA9") = Y.l(z,t);
report(t, z, "I", "SA9") = I.l(z,t);
report(t, z, "a", "SA9") = a.l(z,t);
report(t, z, "C", "SA9") = C.l(t);
report(t, z, "alpha", "SA9") = alpha(z);
report(t, z, "m", "SA9") = m.l(t);
report(t, z, "ko(z)", "SA9") = k0(z);
report(t, z, "v0(z)", "SA9") = v0(z);
report(t, z, "lz(z,t)", "SA9") = lz.l(z,t);

reportt(t, "rho", "SA9") = rho;
reportt(t, "qref", "SA9") = qref(t);
reportt(t, "l0", "SA9") = rho;
reportt(t, "lt", "SA9") = l(t);
reportt(t, "zeta", "SA9") = zeta;

cr("space") = init;

*SA10
init = cr("space");
cr("space") = SA10;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA10") = K.l(z,t);
report(t, z, "R", "SA10") = R.l(z,t);
report(t, z, "V", "SA10") = V.l(z,t);
report(t, z, "Y", "SA10") = Y.l(z,t);
report(t, z, "I", "SA10") = I.l(z,t);
report(t, z, "a", "SA10") = a.l(z,t);
report(t, z, "C", "SA10") = C.l(t);
report(t, z, "alpha", "SA10") = alpha(z);
report(t, z, "m", "SA10") = m.l(t);
report(t, z, "ko(z)", "SA10") = k0(z);
report(t, z, "v0(z)", "SA10") = v0(z);
report(t, z, "lz(z,t)", "SA10") = lz.l(z,t);

reportt(t, "rho", "SA10") = rho;
reportt(t, "qref", "SA10") = qref(t);
reportt(t, "l0", "SA10") = rho;
reportt(t, "lt", "SA10") = l(t);
reportt(t, "zeta", "SA10") = zeta;

cr("space") = init;


*SA11
init = cr("earth");
cr("earth") = SA11;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA11") = K.l(z,t);
report(t, z, "R", "SA11") = R.l(z,t);
report(t, z, "V", "SA11") = V.l(z,t);
report(t, z, "Y", "SA11") = Y.l(z,t);
report(t, z, "I", "SA11") = I.l(z,t);
report(t, z, "a", "SA11") = a.l(z,t);
report(t, z, "C", "SA11") = C.l(t);
report(t, z, "alpha", "SA11") = alpha(z);
report(t, z, "m", "SA11") = m.l(t);
report(t, z, "ko(z)", "SA11") = k0(z);
report(t, z, "v0(z)", "SA11") = v0(z);
report(t, z, "lz(z,t)", "SA11") = lz.l(z,t);

reportt(t, "rho", "SA11") = rho;
reportt(t, "qref", "SA11") = qref(t);
reportt(t, "l0", "SA11") = rho;
reportt(t, "lt", "SA11") = l(t);
reportt(t, "zeta", "SA11") = zeta;

cr("earth") = init;


*SA12
init = cr("earth");
cr("earth") = SA12;

solve noterm maximizing W using NLP;

report(t, z, "K", "SA12") = K.l(z,t);
report(t, z, "R", "SA12") = R.l(z,t);
report(t, z, "V", "SA12") = V.l(z,t);
report(t, z, "Y", "SA12") = Y.l(z,t);
report(t, z, "I", "SA12") = I.l(z,t);
report(t, z, "a", "SA12") = a.l(z,t);
report(t, z, "C", "SA12") = C.l(t);
report(t, z, "alpha", "SA12") = alpha(z);
report(t, z, "m", "SA12") = m.l(t);
report(t, z, "ko(z)", "SA12") = k0(z);
report(t, z, "v0(z)", "SA12") = v0(z);
report(t, z, "lz(z,t)", "SA12") = lz.l(z,t);

reportt(t, "rho", "SA12") = rho;
reportt(t, "qref", "SA12") = qref(t);
reportt(t, "l0", "SA12") = rho;
reportt(t, "lt", "SA12") = l(t);
reportt(t, "zeta", "SA12") = zeta;

cr("earth") = init;

execute_unload "baseramsey.gdx" report, reportt;