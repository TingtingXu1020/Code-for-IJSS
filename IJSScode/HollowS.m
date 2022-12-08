function [NMicporosity, NGmic_compliance_vp, NGmic_compliance_e, Nvoid_r,  Neps_0]...
    = HollowS(num, Micsig, dMicsig, void_r,  Q, Timeinterval, K, mu, Cell_R, scale_up, P_l, Micporosity, nu, R, T)
% Using dMicsig in last step to update the stress as initial guess
NMicsig = Micsig(:,:,num) + dMicsig(:,:,num)*Timeinterval; 
Gmic_compliance_e = HScompliance(K, mu, Micporosity(num), nu);
%% ====================================== Chemical Response ========================================
% Healing response 
% The updated local stress LMicsig enters into healing mechanism directly
LMicsig = Transform(NMicsig, Q(:, :, num));
[dvoid_r, Lmic_compliance_c]  = Healing(LMicsig, void_r(num), Cell_R, scale_up, R, T, P_l);
% Update state variable and Coordinate transform
Gmic_compliance_c = Transform(Lmic_compliance_c, transpose(Q(:, :, num)));
% !!!!!! Approximation of chemical compliance (Penalty method for shear stiffness)
if norm(tensor2matrix(Gmic_compliance_c)) > 0 && rcond(tensor2matrix(Gmic_compliance_c)) < 1E-16
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Gmic_compliance_c(i, j, k, l)= Gmic_compliance_c(i, j, k, l)...
                        + 10^(-6)*norm(tensor2matrix(Gmic_compliance_c))*(KronD(i,k)*KronD(j,l) + KronD(i,l)*KronD(j,k));
                end
            end
        end
    end
end
Nvoid_r = void_r(num) +  dvoid_r *Timeinterval;
if Nvoid_r < 0
    Nvoid_r  = 0;
end
NMicporosity = Nvoid_r^3/Cell_R^3;
NGmic_compliance_e = HScompliance(K, mu, NMicporosity, nu);
NGmic_compliance_vp = Gmic_compliance_c/(1-NMicporosity)...
    + (NGmic_compliance_e - Gmic_compliance_e)/Timeinterval;
Neps_0 = doubledotft(Gmic_compliance_c/(1-NMicporosity), P_l*eye(3))*NMicporosity;
end