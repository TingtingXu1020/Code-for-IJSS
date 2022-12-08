function [SSumMB_evp, SSumMPhie, SSumMB_e, SSumB, SSumPhi] =...
    Acco(Gmic_compliance_e,  Gmic_compliance_evp,...
    num, Mac_compliance_evp, Interaction_e, Interaction_evp, weight, ...
    dMiceps_0, Micsig, Macsig, dE0, Timeinterval)
B =  doubledotff(inversegeneral(Gmic_compliance_evp(:, :, :, :, num)  +...
    Interaction_evp), (Mac_compliance_evp + Interaction_evp));
de0 = -1/Timeinterval * doubledotft(Gmic_compliance_e(:, :, :, :, num), Micsig(:, :, num)) + ...
    dMiceps_0(:, :, num);
temp1 = dE0 - de0 -1/Timeinterval * doubledotft(Interaction_e, (Macsig - Micsig(:, :, num)));
Phi = doubledotft(inversegeneral(Gmic_compliance_evp(:, :, :, :, num)  +...
    Interaction_evp), temp1);

if isnan(B)
    error('nan')
end
SSumMB_evp  =  doubledotff(Gmic_compliance_evp(:, :, :, :, num), B) * weight(num);
SSumMPhie  =  (doubledotft(Gmic_compliance_evp(:, :, :, :, num), Phi) + de0) * weight(num);
SSumMB_e  =  doubledotff(Gmic_compliance_e(:, :, :, :, num), B) * weight(num);
SSumB =  B * weight(num);
SSumPhi = Phi * weight(num);
end