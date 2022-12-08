function [ tauini_h,  dvoid_r, Gmic_compliance_sl, eps_back ] = ...
    Sliding(Micsig, aglobal,  tauini, approx, n_slip, R, T)
% Parameters                                                  
slipsys = 6;  
Q =  8e6; % Activation energy  mJ/mol
gamma_0 = 1*exp(-Q/R/T); 
% Initial set 
% All under global coordinate system
dgamma_slip = 0;
d_sumgamma = zeros(3,3);
Gmic_compliance_sl = zeros(3, 3, 3, 3);
for slip = 1 : slipsys
    tauslip = doubledottt(Micsig, aglobal(:, :, slip));
    Gmic_compliance_sl = Gmic_compliance_sl + gamma_0*(abs(tauslip/tauini))^(n_slip - 1)/ tauini ...
        * outerproducttt(aglobal(: ,:, slip), aglobal(: ,:, slip));
    dgamma_slip = dgamma_slip + gamma_0* sign(tauslip)*(abs(tauslip / tauini)^(n_slip));
    d_sumgamma = d_sumgamma + gamma_0* sign(tauslip)*(abs(tauslip / tauini)^(n_slip)) * aglobal(: ,:, slip);
end
switch approx
    case 1 % tangent scheme
        eps_back = (1 - n_slip)*d_sumgamma;
        Gmic_compliance_sl = n_slip * Gmic_compliance_sl;
    case 2 % secant scheme
        eps_back = zeros(3, 3);
        Gmic_compliance_sl = 1* Gmic_compliance_sl;
end
tauini_h = tauini;
% Intact elastic response
dvoid_r = 0;
end