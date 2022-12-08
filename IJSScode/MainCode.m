delete(gcp('nocreate'));
num_core = 8;
c = parcluster('local');
c.NumWorkers = num_core;
p = parpool(c);
clc
close all
clear
dbstop if error

%===================================================================%
%                        Inclusion generation                            %
%===================================================================%
hydro = 3;
porosity0 = 0.4;
Cell_R = 0.200;
%% case 1 (dry sample -- only have single grain inclusion and pore)
switch hydro
    case 1 % (dry sample -- only have single grain inclusion and pore)
        IntpointsA = [0; 0; 0];
        IntpointsB = [500; 500; 500];
        f_hs0 = 0.0;
        f_p0 = porosity0;
        f_s0 = 1.0 - f_p0;
        
    case 2 % (saturated sample only considering pressure solution mechanism)
        IntpointsA = [500; 500; 500];
        IntpointsB = [0; 0; 0];
        f_hs0 = 1.0;
        f_s0 = 0;
        f_p0 = 0;
        % Shape parameters
        a_0 = (Cell_R^3*porosity0)^(1/3); % Average pore radius
        rVar = 1e-5; % variance of pore radius
        
    case 3 % (saturated sample considering both mechanisms with different degree of saturation)
        IntpointsA = [500; 500; 500];
        IntpointsB = IntpointsA;
        DoS = 1; % control the dry condition
        f_p0 = (1-DoS)*porosity0;
        f_compo = 1 - porosity0;
        f_s0 = f_compo/2;
        f_hs0 =  1 - f_p0 - f_s0;     
        % Shape parameters
        a_0 = (Cell_R^3*(DoS * porosity0/f_hs0))^(1/3); % Average pore radius
        rVar = 1e-5; % variance of pore radius
end
cellnumberA = IntpointsA(1);
cellnumberB = IntpointsB(1);
cellnumber = cellnumberA + cellnumberB;
method = 1; % 1 indicates Jiang's method, while 2 indicates Pouya's method
weight = zeros(cellnumber, 1);
if norm(IntpointsA) > 0
    [angle_A, rho_A, R_A, number_A] = Generation(IntpointsA(1), IntpointsA(2), IntpointsA(3), method);
    cellnumberA = number_A;
    weight(1:cellnumberA) =  rho_A * f_hs0;
    Q(:, :, 1 : cellnumberA) = R_A;
end
if norm(IntpointsB) > 0
    [angle_B, rho_B, R_B, number_B] = Generation(IntpointsB(1), IntpointsB(2), IntpointsB(3), method);
    cellnumberB = number_B;
    weight((cellnumberA + 1) : cellnumber) = rho_B * f_s0;
    Q(:, :, (cellnumberA + 1) : cellnumber) = R_B;
end

% Set the initial value
f_hs = f_hs0;
f_s = f_s0;
f_p = f_p0;

%===================================================================%
%                        Input parameters                            %
%===================================================================%
% Elastic parameter
E = 2.5e4; % Young's modulus (MPa)
nu = 0.37;  % Possion's ratio
K = E/3/(1 - 2*nu);
lambda = E*nu /(1 + nu)/(1 - 2*nu);
mu = E/2/(1+nu);

mu_sc = 3*mu*(1-2*porosity0)/(3-porosity0);
K_sc = 4*mu*(1-2*porosity0)*(1-porosity0)/porosity0/(3-porosity0);
E_sc = 9*K_sc*mu_sc/(3*K_sc + mu_sc);
nu_sc = (3*K_sc - 2* mu_sc )/2/(3*K_sc + mu_sc);
% Damage parameters
sigtensile = 200;
% Sliding parameter
n_slip   = 8.6;
tauini  = 5;
n_eff = 1.0; % effective tunning variable
approx = 1; % 1 indicates tangent approximation, while 2 indicates secant approximation
% Healing parameter (in order to scale the influence of healing)
scale_up = 1e0; % in order to scale the influence of healing
% Temperature
R = 8.3145e3;   % mJ/mol/K
T = 273 + 22;    % K
% Numerical parameter
int_e = [16; 16];
int_evp = [16; 100];

%===================================================================%
%                        Boundary Conditions                        %
%===================================================================%
% Boundary conditions
inc = 1;
Timeinterval = 10; % Dt
sig_c = 0;    % Confining stress
control_matrix = [  -1 0 0; 0 0 0; 0 0 0];
Lstep = 24*3600*4; % Number of total loading steps
P_l = 0.1; % Pressure of liquid (1atm)
P_g = 0.1; % Pressure of gas (1atm)
control = 2;

switch control
    case 1 % Creep Test (Uniaxial condition)
        sig_a = -16.5; % Applied Axial stress (MPa)
        Macsig = eye(3) * sig_c;
        [row, col] = find(control_matrix);
        Macsig(row, col) = sig_a;
        Maceps = zeros(3, 3);
        dMaceps = zeros(3, 3);
        dMacepse = dMaceps;
        dMacepsvp = zeros(3, 3);
    case 2 % Creep Test (Oedemetric condition)
        sig_a = -2.1; % Applied Axial stress (MPa)
        Macsig = eye(3) * sig_c;
        [row, col] = find(control_matrix);
        Macsig(row, col) = sig_a;
        Spiers = xlsread('Spiers1990');
        Maceps = [-Spiers(1,6)/100, 0, 0; 0, 0, 0; 0, 0, 0] ;
        dMaceps = zeros(3, 3);
        dMacepse = dMaceps;
        dMacepsvp = zeros(3, 3);
    case 3 % Mixed Control (Strain control Uniaxial test)
        dMaceps11 = 5e-7; % s^(-1) strain rate
        dMaceps = dMaceps11 * control_matrix;
        Maceps = zeros(3, 3);
        Macsig = eye(3) * sig_c;
        dMacepse = dMaceps;
        dMacepsvp = zeros(3, 3);
end

%===================================================================%
%                     Initialization                         %
%===================================================================%
% Microscopic state variables
dMicepse(1:3,1:3, 1:cellnumber) = 0;
dMicepsvp(1:3,1:3, 1:cellnumber) = 0;
Micepse(1:3,1:3, 1:cellnumber) = 0;
Micepsvp(1:3,1:3, 1:cellnumber) = 0;
dMiceps_0 (1:3, 1:3, 1:cellnumber) = 0;
Alldstrain_vpA(:, :) = zeros(3, 3);
Alldstrain_vpB(:, :) = zeros(3, 3);
Alldstrain_vpC(:, :) = zeros(3, 3);

% Void allocation
void_r(1:cellnumber) = 0;
connect_porosity(inc) = 0;
if cellnumberA > 0
    void_r(1:cellnumberA) = lognrnd(log(a_0^2/sqrt(rVar + a_0^2)), sqrt(log(rVar/a_0^2 + 1)), cellnumberA, 1);
    connect_porosity(inc) = sum(void_r.^3*weight)/Cell_R^3;
end
nonconnect_porosity(inc) = f_p; % the fraction of the volume of voids over the total volume
Macporosity(inc) = connect_porosity(inc) + nonconnect_porosity(inc);

% Macroscopic state variables
dMaceps_0 = zeros(3, 3);
Maceps_v(inc) = trace(Maceps);

% Storage matrix of boundary conditions and time
Macsig_save(:, :, inc) = Macsig;
Maceps_save(:, :, inc) = Maceps;
T_save(inc) = 0;
%% First step
[Jd, Js, Ja, Jm] = Identity();

% Microscopic compliance for the intact material
Micporosity(:) = void_r.^3/Cell_R^3;
Gmic_compliance_vp(1:3, 1:3, 1:3, 1:3, 1:cellnumber) = 0;
Gmic_compliance_e(1:3, 1:3, 1:3, 1:3, 1:cellnumber) = 0;
for num =1 : cellnumber
    Gmic_compliance_vp(:, :, :, :, num)  = zeros(3, 3, 3, 3);
    Gmic_compliance_e(:, :, :, :, num)  = Intactcompliance(E, nu);
end

% Initial Macroscopic compliance (Intact material)
Mac_compliance_e = Intactcompliance(E_sc, nu_sc);
Mac_compliance_vp = zeros(3, 3, 3, 3);

switch control
    case 1
        dMacsig = zeros(3, 3);
    case 2
        dMacsig = zeros(3, 3);
    case 3
        dMacsig = doubledotft(inversegeneral(Mac_compliance_e), dMaceps);
end

dMicsig(1:3, 1:3, 1:cellnumber) = 0;
Micsig(1:3, 1:3, 1:cellnumber) = 0;
dMiceps(1:3, 1:3, 1:cellnumber) = 0;
Miceps(1:3, 1:3, 1:cellnumber) = 0;
for num = 1 : cellnumber
    dMicsig(:, :, num) = dMacsig;
    Micsig(:, :, num) = Macsig;
    dMiceps(:, :, num) = dMaceps;
    Miceps(:, :, num) = Maceps;
end

% Dislocation calculation
%%% == sliding and normal vector == %%%
A1 = zeros(3, 3); A2 = zeros(3, 3); A3 = zeros(3, 3);
N1 = [0 1 1]' /sqrt(2);
N2 = [1 0 1]' / sqrt(2);
N3 = [-1 -1 0]' / sqrt(2);
N4 = [0 -1 1]' / sqrt(2);
N5 = [-1 0 1]' / sqrt(2);
N6 = [-1 1 0]' / sqrt(2);
M1 = -N4; M2 = -N5; M3 = -N6;
M4 = -N1; M5 = -N2; M6 = -N3;
% Schmid tensor
for i = 1:3
    for j = 1:3
        A1(i, j) = 0.5 * (M1(i) * N1(j) + M1(j) * N1(i));
        A2(i, j) = 0.5 * (M2(i) * N2(j) + M2(j) * N2(i));
        A3(i, j) = 0.5 * (M3(i) * N3(j) + M3(j) * N3(i));
    end
end
A_slide(:, :, 1) = A1; A_slide(:, :, 2) = A2; A_slide(:, :, 3) = A3;
A_slide(:, :, 4) = A1; A_slide(:, :, 5) = A2; A_slide(:, :, 6) = A3;

% Iteraction Check
Titer = 1;
Res_sig = zeros(cellnumber,1);

%% Main code
RepI = [];
Hollownum = [];
Singlenum = [];
if cellnumberA > 0
    Hollownum = 1 : cellnumberA;
end
if cellnumberB > 0
    Singlenum = cellnumberA + 1 : cellnumber;
end
tnum = [Hollownum, Singlenum];

while inc < Lstep % Total loading step
    T_save(inc) = Timeinterval*inc;
    disp(['==========',num2str(inc),'==========='])
    % save matrix
    Macsig_save(:, :, inc) = Macsig;
    Maceps_save(:, :, inc) = Maceps;
    Mac_compliance_esave(:,:,:,:,inc) = Mac_compliance_e;
    Mac_compliance_vpsave(:,:,:,:,inc) = Mac_compliance_vp;
    %% Update Macsig and Maceps and loading step
    inc = inc + 1;
    Titer = 1;
    if inc == 2
        Tnum = 5;
    else
        Tnum = 50;
    end
    %% Outermost loop for Micsig Iteration
    while Titer < Tnum
        %===================================================================%
        %                        Inclusion behavior  (All state variables are from last step)                      %
        %===================================================================%
        %% ==================================== Hollow sphere inclusion Response ===============================
        if ~isempty(Hollownum)
            parfor inum = 1 : length(Hollownum)
                num = tnum(inum);
                [NMicporosity, NGmic_compliance_vp, NGmic_compliance_e, Nvoid_r,  Neps_0]...
                    = HollowS(num, Micsig, dMicsig, void_r, Q, Timeinterval, K, mu, Cell_R,...
                    scale_up, P_l, Micporosity, nu, R, T);
                ovoid_r(inum) = Nvoid_r;
                oMicporosity(inum) = NMicporosity;
                oGmic_compliance_vp(:, :, :, :, inum) = NGmic_compliance_vp;
                oGmic_compliance_e(:, :, :, :, inum) = NGmic_compliance_e;
                oGmic_compliance_evp(:, :, :, :, inum) = 1/Timeinterval * NGmic_compliance_e + NGmic_compliance_vp;
                oeps_0(:, :, inum) = Neps_0;
            end
        end
        %% ====================================== Single-crystal response =====================================
        if ~isempty(Singlenum)
            parfor inum = length(Hollownum) + 1 : length(tnum)
                num = tnum(inum);
                [Nvoid_r, Ntauini,  NGmic_compliance_vp, Neps_0] ...
                    = SingleI(num, Micsig, dMicsig, Timeinterval, A_slide, Q,  tauini, void_r,  approx, n_slip, R, T);
                ovoid_r(inum) = Nvoid_r;
                oGmic_compliance_vp(:, :, :, :, inum) = NGmic_compliance_vp;
                oGmic_compliance_e(:, :, :, :, inum)  =  Intactcompliance(E, nu);
                oGmic_compliance_evp(:, :, :, :, inum) = 1/Timeinterval * Intactcompliance(E, nu) + NGmic_compliance_vp;
                oeps_0(:, :, inum) = Neps_0;
                oMicporosity(inum) = 0;
            end
        end
        
        void_r = ovoid_r;
        Micporosity = oMicporosity;
        dMiceps_0 = oeps_0;
        Gmic_compliance_vp = oGmic_compliance_vp;
        Gmic_compliance_e = oGmic_compliance_e;
        Gmic_compliance_evp = oGmic_compliance_evp;
        
        % Average for the first step
        if inc == 2
            Compliance_vp = zeros(3,3,3,3);
            Compliance_e = zeros(3,3,3,3);
            for num = 1 : cellnumber
                Compliance_vp = Compliance_vp + Gmic_compliance_vp(:, :, :, :, num) * weight(num);
                Compliance_e = Compliance_e + Gmic_compliance_e(:, :, :, :, num) * weight(num)  ;
                dMaceps_0 = dMaceps_0 + dMiceps_0(:, :, num) * weight(num);
            end
            Mac_compliance_vp = Symmetry(Compliance_vp);
            Mac_compliance_e = Symmetry(Compliance_e);
            dE0 = -1/Timeinterval * doubledotft(Mac_compliance_e, Macsig) + dMaceps_0;
        end
        %===================================================================%
        %                        Homogenization to RVE Scale                         %
        %===================================================================%
        %% innermost loop for Macroscopic compliances and Eshelby tensors
        j_iter = 300;
        Mac_compliance_evp = Mac_compliance_vp + Mac_compliance_e/Timeinterval;
        [Mac_compliance_evp, Mac_compliance_e, dE0, Interaction_evp, Interaction_e,...
            Error, Div, B, Phi] = Maccompliance( ...
            Gmic_compliance_e, Gmic_compliance_evp, Js,   j_iter,...
            Mac_compliance_e, Mac_compliance_vp, Mac_compliance_evp,  weight, f_p , Cell_R,...
            Macsig, int_evp, tnum, Timeinterval,  dMiceps_0, Micsig, dE0, n_eff, porosity0);
        Mac_compliance_vp = Mac_compliance_evp - 1/Timeinterval*Mac_compliance_e;
        dMaceps_0 = dE0 + 1/Timeinterval * doubledotft(Mac_compliance_e, Macsig); 
        if Div == 1
            disp('diverse')
        end
        
        %% Upscale Macstrain or Macstress
        M_tensor = Mac_compliance_evp;
        N_matrix =  dMaceps_0;
        % Sig = M^(-1) * (dE -N);
        switch control
            case 1 % uniaxial creep test
                dMacsig = zeros(3,3);
                dMaceps = doubledotft(Mac_compliance_evp, Macsig) + dMaceps_evp0;
            case 2 % oedemetric creep test
                M_matrix = tensor2matrix(M_tensor);
                N_vector = tensor2vector(N_matrix);
                ReduM = M_matrix(2:6, 2:6);
                ReduC = inv(ReduM);
                tempLeft = zeros(5, 1);
                for i = 1 : 5
                    tempLeft(i) = -M_matrix(i+1, 1)*Macsig(1,1);
                end
                tempstress = ReduC*(tempLeft - N_vector(2 : 6));
                tempstress = [Macsig(1,1);tempstress];
                Macsig = vector2tensor(tempstress);
                dMaceps = doubledotft(Mac_compliance_evp, Macsig) + dMaceps_0;
            case 3 % uniaxial loading
                M_matrix = tensor2matrix(M_tensor);
                N_vector = tensor2vector(N_matrix);
                dMaceps(1, 1) = dMaceps11 * control_matrix(1,1);
                dMacsig(1, 1) = 1/M_matrix(1,1) * (dMaceps(1, 1) - N_vector(1));
                dMacsig(2, 2) = 0;
                dMacsig(3, 3) = 0;
                dMaceps = doubledotft(M_tensor, dMacsig) + N_matrix;
                dMaceps(1, 1) = dMaceps11 * control_matrix(1,1);
        end
        %% Update local stress and stress rate
        if ~isempty(RepI)
            Micsig(1:3, 1:3, RepI) = 0;
            dMicsig(1:3, 1:3, RepI) = 0;
        end
        Macsig = Macsig + dMacsig*Timeinterval;
        Maceps = Maceps + dMaceps*Timeinterval;
        avesig = zeros(3, 3);
        for inum = 1: length(tnum)
            num = tnum(inum);
            NMicsig = doubledotft(B(:, :, :, :, num), Macsig) + Phi(:, :, num);
            Res_sig(num) = norm(NMicsig - Micsig(:, :, num))/ norm(Micsig(:, :, num));
            dMicsig(:, :, num) = (NMicsig - Micsig(:, :, num))/Timeinterval;
            Micsig(:, :, num) = NMicsig;
            avesig = avesig + Micsig(:, :, num) * weight(num);
        end
        errorsig = avesig - Macsig;
        %% check updated strain
        for inum = 1: length(tnum)
            num = tnum(inum);
            dMicepse(:, :, num) = doubledotft(Gmic_compliance_e(:, :, :, :, num), dMicsig(:, :, num) );
            dMicepsvp(:, :, num) = doubledotft(Gmic_compliance_vp(:, :, :, :, num), Micsig(:, :, num));
            dMiceps(:, :, num) = dMicepse(:, :, num) + dMicepsvp(:, :, num) + dMiceps_0(:, :, num);
            Miceps(:, :, num) =  dMiceps(:, :, num)*Timeinterval + Miceps(:, :, num);
        end
        %% Convergence Check
        if mean(Res_sig) < 0.06
            break
        end
        Titer = Titer + 1;
        %         if Titer == Tnum
        %             disp('Current step is diverse');
        %         end
        mean(Res_sig)
    end
    %% End iteration and update
    connect_porosity(inc) = sum(void_r.^3 * weight)/Cell_R^3; %% assume the weight only transfer between
    % single grain and pore, do not influence the initial volume fraction
    % of pore
    if connect_porosity(inc) > connect_porosity(inc-1)
        error('error')
    end
    %% Replacement occurs and use the connect porosity together with non-connect porosity to update the total porosity
    temp = zeros(cellnumber, 1);
    for num = cellnumberA + 1 : cellnumber
        temp(num) = max(eig(Micsig(:, :, num)));
    end
    RepI = [RepI; find(temp > sigtensile)];
    tnum = sort(setdiff(tnum, RepI));
    Singlenum = sort(setdiff(Singlenum, RepI));
    
    if ~isempty(RepI)
        f_p = f_p + sum(weight(RepI));
    end
    
    weight(RepI) = 0;
    
    if abs(sum(weight) + f_p - 1.0 ) > 1e-6
        error('fraction error')
    end
    %
    if f_p > 1/2
        error('too many broken inclusions')
    end
        
    nonconnect_porosity(inc) = f_p;
    Maceps_v(inc) = Maceps_v(inc-1) +  trace(dMaceps)*Timeinterval;
    Macporosity(inc) = Macporosity(inc-1) + trace(dMaceps)*(1-Macporosity(inc-1))*Timeinterval;
end
if DoS < 1E-2
DoS = 0;
end
Active = hydro-1; % 0 represents only pressure solution, 1 represents only dislocation 
save(['DoS',num2str(DoS) , 'Time',num2str(T_save(end)/24/3600), 'Siga', num2str(sig_a), 'Micsig.mat'],'Micsig');
save(['DoS',num2str(DoS) ,'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'T.mat'],'T_save');
save(['DoS',num2str(DoS) , 'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'Macporosity.mat'],'Macporosity');
save(['DoS',num2str(DoS) ,'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'Maceps_v.mat'],'Maceps_v');
save(['DoS',num2str(DoS) ,'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'Maceps_save.mat'],'Maceps_save');
save(['DoS',num2str(DoS) ,'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'Macsig_save.mat'],'Macsig_save');
save(['DoS',num2str(DoS) ,'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'Mac_compliance_esave.mat'],'Mac_compliance_esave');
save(['DoS',num2str(DoS) ,'Time',num2str(T_save(end)/24/3600),  'Siga', num2str(sig_a), 'Mac_compliance_vpsave.mat'],'Mac_compliance_vpsave');

Macporosity(1) = Macporosity(2);
Maceps_save(:, :, 1) = Maceps_save(:, :, 2);
Plot(cellnumber, Micsig, inc, Maceps_save,...
    Macsig_save, T_save,  Macporosity, Maceps_v,... 
    nonconnect_porosity, cellnumberA, connect_porosity);

