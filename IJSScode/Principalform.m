% Principal direction derivation
syms D1 D2 D3 sigtensile E eta presig preeps11 preeps22 preeps33 preeps12 preeps13 preeps23
syms  lambda mu eps11 eps22 eps33 eps12 eps13 eps23 Timeinterval

C11 = (lambda + 2*mu) * (1-D1)^(1/2) * (1-D1)^(1/2);
C12 = lambda * (1-D1)^(1/2) * (1-D2)^(1/2);
C13 = lambda * (1-D1)^(1/2) * (1-D3)^(1/2);
C22 = (lambda + 2*mu) *(1-D2)^(1/2) * (1-D2)^(1/2);
C33 = (lambda + 2*mu) *(1-D3)^(1/2) * (1-D3)^(1/2);
C23 = lambda * (1-D2)^(1/2) * (1-D3)^(1/2);
C44 = 2 * mu * (1-D2)^(1/2) * (1-D3)^(1/2);
C55 = 2 * mu * (1-D1)^(1/2) * (1-D3)^(1/2);
C66 = 2 * mu * (1-D1)^(1/2) * (1-D2)^(1/2);
C = [C11 C12 C13 0 0 0
    C12 C22 C23 0 0 0
    C13 C23 C33 0 0 0
    0 0 0 C44 0 0
    0 0 0 0 C55 0
    0 0 0 0 0 C66];
eps = [eps11; eps22; eps33; eps23; eps13; eps12];

% syms sig11  sig22  sig33  sig12  sig13  sig23
% stress = [sig11 sig12 sig13
%     sig12 sig22 sig23
%     sig13 sig23 sig33];

sig = C*eps;
stress = [sig(1) sig(6) sig(5);
    sig(6) sig(2) sig(4)
    sig(5) sig(4) sig(3)];
strain = [eps(1) eps(6) eps(5);
    eps(6) eps(2) eps(4)
    eps(5) eps(4) eps(3)]; 
lmicdamage = [D1; D2; D3];
e = sym('e', [3,3]);
e = real(e);
Fd = e' *stress*e - ((sigtensile -  E *e' *strain*e + eta*sigtensile) /eta);


Jd11 = diff(Fd(1), D1);
Jd12 = diff(Fd(1), D2);
Jd13 = diff(Fd(1), D3);
Jd21 = diff(Fd(2), D1);
Jd22 = diff(Fd(2), D2);
Jd23 = diff(Fd(2), D3);
Jd31 = diff(Fd(3), D1);
Jd32 = diff(Fd(3), D2);
Jd33 = diff(Fd(3), D3);


Jeps11 = diff(Fd(1), eps11);
Jeps12 = diff(Fd(1), eps22);
Jeps13 = diff(Fd(1), eps33);
Jeps21 = diff(Fd(2), eps11);
Jeps22 = diff(Fd(2), eps22);
Jeps23 = diff(Fd(2), eps33);
Jeps31 = diff(Fd(3), eps11);
Jeps32 = diff(Fd(3), eps22);
Jeps33 = diff(Fd(3), eps33);