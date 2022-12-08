syms D1 D2 D3 sigtensile E eta presig preeps11 preeps22 preeps33 preeps12 preeps13 preeps23
syms  lambda mu eps11 eps22 eps33 eps12 eps13 eps23 Timeinterval

ncosine = eye(3);
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

syms sig11 sig22 sig33 sig23 sig13 sig12
stress = [sig11 sig12 sig13;
    sig21 sig22 sig23
    sig31 sig32 sig33];
strain = [eps(1) eps(6) eps(5);
    eps(6) eps(2) eps(4)
    eps(5) eps(4) eps(3)]; 
Fd1 = ncosine(1, :) * stress * transpose(ncosine(1, :)) - ...
    (( (sigtensile -  E *(ncosine(1, :) * strain * transpose(ncosine(1, :))) + eta*sigtensile) /eta));
Fd2 = ncosine(2, :) * stress * transpose(ncosine(2, :)) - ...
    (( (sigtensile -  E *(ncosine(2, :) * strain * transpose(ncosine(2, :))) + eta*sigtensile) /eta));
Fd3 = ncosine(3, :) * stress * transpose(ncosine(3, :)) - ...
    (( (sigtensile -  E *(ncosine(3, :) * strain * transpose(ncosine(3, :))) + eta*sigtensile) /eta));

