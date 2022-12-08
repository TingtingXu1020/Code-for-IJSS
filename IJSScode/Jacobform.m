
syms  lambda mu Timeinterval lmicdamage1 lmicdamage2 lmicdamage3

lmicdamage = [lmicdamage1; lmicdamage2; lmicdamage3];

C11 = (lambda + 2*mu) * (1- lmicdamage1)^(1/2) * (1- lmicdamage1)^(1/2);
C12 = lambda * (1- lmicdamage1)^(1/2) * (1- lmicdamage2)^(1/2);
C13 = lambda * (1- lmicdamage1)^(1/2) * (1- lmicdamage3)^(1/2);
C22 = (lambda + 2*mu) *(1- lmicdamage2)^(1/2) * (1- lmicdamage2)^(1/2);
C33 = (lambda + 2*mu) *(1- lmicdamage3)^(1/2) * (1- lmicdamage3)^(1/2);
C23 = lambda * (1- lmicdamage2)^(1/2) * (1- lmicdamage3)^(1/2);
C44 = 2 * mu * (1- lmicdamage2)^(1/2) * (1- lmicdamage3)^(1/2);
C55 = 2 * mu * (1- lmicdamage1)^(1/2) * (1- lmicdamage3)^(1/2);
C66 = 2 * mu * (1- lmicdamage1)^(1/2) * (1- lmicdamage2)^(1/2);
C = [C11 C12 C13 0 0 0
    C12 C22 C23 0 0 0
    C13 C23 C33 0 0 0
    0 0 0 C44 0 0
    0 0 0 0 C55 0
    0 0 0 0 0 C66];

syms leps11 leps22 leps33 leps23 leps13 leps12
syms dleps11 dleps22 dleps33 dleps23 dleps13 dleps12
leps = [leps11;leps22;leps33;leps23;leps13;leps12];
dleps = [dleps11;dleps22;dleps33;dleps23;dleps13;dleps12];
syms presig11  presig22  presig33  presig12  presig13  presig23
 presig = [presig11; presig22; presig33; presig23; presig13; presig12];
% stress = [sig11 sig12 sig13
%     sig12 sig22 sig23
%     sig13 sig23 sig33];

strain = [eps(1) eps(6) eps(5);
   leps(6) eps(2) eps(4)
   leps(5) eps(4) eps(3)];
sig = presig + C*dleps;
stress = [sig(1) sig(6) sig(5);
    sig(6) sig(2) sig(4)
    sig(5) sig(4) sig(3)];
syms sigtensile eta E d_c
xi = lmicdamage/d_c;
dhardening = (1 + eta * xi./(1 + xi.^2));
Fd =  diag(stress) - sigtensile*dhardening;
Jd = jacobian(Fd, lmicdamage);


