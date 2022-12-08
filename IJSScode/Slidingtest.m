%function [ tauini_h, Lmic_compliance_e, dvoid_r, Gmic_compliance_sl ] = Sliding(Micsig, aglobal,  Mic_compliance_intact, tauini)
% Elastic parameter (Ref. Liang 2006)
E = 1E4; % Young's modulus
nu = 0.31;  % Possion's ratio
K = E/3/(1 - 2*nu);
lambda = E*nu /(1 + nu)/(1 - 2*nu);
mu = E/2/(1+nu);
Timeinterval = 1;
% Sliding parameter
tauini   = 5;
% Intact elastic stiffness
Lmic_stiffness_e = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l=1:3
                Lmic_stiffness_e(i,j,k,l) = Lmic_stiffness_e(i,j,k,l)...
                    + lambda *KronD(i,j)*KronD(k,l) + mu*(KronD(i,k)*KronD(j,l) + KronD(i,l)*KronD(j,k));
            end
        end
    end
end
Mic_compliance_intact = inversegeneral(Lmic_stiffness_e);
% Initial condition
Q = [ -0.6103   -0.7918    0.0234
    0.7862   -0.6091   -0.1048
    0.0973   -0.0456    0.9942];
Micsig = [5, 0, 0; 0, 0, 0; 0, 0, 0];
LMicsig = zeros(3, 3);
dMicsig = [0, 0, 0; 0, 0, 0; 0, 0, 0]*1E-5;
Miceps = zeros(3, 3);
% Parameters
gamma_0 = 1;
n_slip        = 4.04;
slipsys = 6;
% Dislocation calculation
%%% == sliding and normal vector == %%%
A1 = zeros(3,3); A2 = zeros(3,3); A3 = zeros(3,3);
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
N(:, :, 1) =  N1*N1'; N(:, :, 2) =  N2*N2'; N(:, :, 3) =  N3*N3'; 
N(:, :, 4) =  N4*N4'; N(:, :, 5) =  N5*N5'; N(:, :, 6) =  N6*N6'; 
A(:, :, 1) = A1; A(:, :, 2) = A2; A(:, :, 3) = A3;
A(:, :, 4) = A1; A(:, :, 5) = A2; A(:, :, 6) = A3;
aglobal(1:3, 1:3, 1:slipsys) = 0; nglobal(1:3, 1:3, 1:slipsys) = 0;
for slip = 1 : slipsys
    aglobal(:, :, slip) = Q*A(:, :, slip)*transpose(Q);
    nglobal(:, :, slip) = Q*N(:, :, slip)*transpose(Q);
end
d_s = [0.1 0.2 0.3 0.2 0.4 0.5];
T = 0;
alpha = 0.1;
Y_0 = 0;
s_d = 1;
S_d = 120;
m_d = 4.05;
while norm(Miceps) < 1
    Micsig = Micsig + dMicsig* Timeinterval;
    % Initial set
    dgamma_slip = 0;
    d_sumgamma = zeros(3,3);
    Gmic_compliance_sl = zeros(3, 3, 3, 3);
    for slip = 1 : slipsys
        tauslip = doubledottt(Micsig, aglobal(:, :, slip));
        lambda_s = gamma_0 * (abs(tauslip / tauini)^(n_slip)); % viscoplastic multiplier
        dgamma_slip = dgamma_slip +  sign(tauslip) * lambda_s/sqrt(1 - d_s(slip));
        d_sumgamma = d_sumgamma + lambda_s/sqrt(1 - d_s(slip)) * (sign(tauslip) *  aglobal(: ,:, slip) ...
            + alpha * nglobal(:, :, slip) *sum(d_s));
        Gmic_compliance_sl = Gmic_compliance_sl + gamma_0* ...
            (abs(tauslip/tauini))^(n_slip - 1)/ tauini /sqrt(1 - d_s(slip)) * (outerproducttt(aglobal(: ,:, slip), aglobal(: ,:, slip)) + ...
            alpha  *sum(d_s) * outerproducttt(nglobal(: ,:, slip), aglobal(: ,:, slip)));
    end
    dMicepssl = doubledotft(Gmic_compliance_sl, Micsig);
    % Intact elastic response
    Lmic_compliance_e = Mic_compliance_intact;
    dvoid_r = 0;
    Miceps = Miceps + dMicepssl + doubledotft(Lmic_compliance_e, dMicsig) * Timeinterval ;
    T = T + 1;
    time(T) = T;
    X(T) = Miceps(1, 1);
    y1(T) = Micsig(1, 1);
    y2(T) = Micsig(2, 2);
    y3(T) = Micsig(3, 3);
end
plot(time, X);