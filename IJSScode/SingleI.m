function [Nvoid_r, Ntauini,  NGmic_compliance_vp, Neps_0] ...
    = SingleI(num, Micsig,dMicsig, Timeinterval, A_slide, Q,  tauini, void_r,  approx, n_slip, R, T)
% Using dMicsig in last step to update the stress as initial guess
NMicsig = Micsig(:,:,num) + dMicsig(:,:,num)*Timeinterval; 
% Slip mechanism
slipsys = 6;
aglobal(1:3, 1:3, 1:slipsys) = 0;
for slip = 1 : slipsys
    aglobal(:, :, slip) = Transform(A_slide(:, :, slip), Q(:, :, num));
end
[ tauini_h,  dvoid_r, Gmic_compliance_sl, eps_back ] =...
    Sliding( NMicsig, aglobal, tauini, approx, n_slip, R, T);
% Update state variable and Coordinate transform
Ntauini = tauini_h;
% !!!!!! Approximation of slide compliance (Penalty method for compressive stiffness)
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Gmic_compliance_sl(i, j, k, l) = Gmic_compliance_sl(i, j, k, l) ...
                    + 10^(-6)*norm(tensor2matrix(Gmic_compliance_sl))*(KronD(i,k)*KronD(j,l) + KronD(i,l)*KronD(j,k));
            end
        end
    end
end
Nvoid_r = void_r(num) +  dvoid_r * Timeinterval;
NGmic_compliance_vp = Gmic_compliance_sl;
Neps_0 =  eps_back;
end