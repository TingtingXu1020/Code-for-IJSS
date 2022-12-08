function         [Mac_compliance_evp, Mac_compliance_e, dE0, Interaction_evp, Interaction_e,...
    Error, Div, B, Phi] = Maccompliance( ...
    Gmic_compliance_e, Gmic_compliance_evp, Js,   j_iter,...
    Mac_compliance_e, Mac_compliance_vp, Mac_compliance_evp,  weight, f_p , Cell_R,...
    Macsig, int_evp, tnum, Timeinterval,  dMiceps_0, Micsig, dE0, n_eff, porosity0)
j = 1;
while j < j_iter
    Mac_compliance_eff = n_eff * Mac_compliance_vp + Mac_compliance_e/Timeinterval;
    Hill_evp = Hill_P(Cell_R*[1;1;1], int_evp(1), int_evp(2), inversegeneral(Mac_compliance_eff));
    if isnan(Hill_evp)
        error('nan')
    end
    S_evp = doubledotff(Symmetry(Hill_evp), inversegeneral(Mac_compliance_evp));
    S_eff = doubledotff(Symmetry(Hill_evp), inversegeneral(Mac_compliance_eff));
    Interaction_evp =  doubledotff(doubledotff(inversegeneral(Js - S_eff), S_eff), Mac_compliance_eff);
    Interaction_e =  doubledotff(doubledotff(inversegeneral(Js - S_evp), S_evp), Mac_compliance_e);
    if isnan(Interaction_evp)
        error('nan')
    end
    parfor inum = 1 : length(tnum)
        num =  tnum(inum);
        [SSumMB_evp, SSumMPhie, SSumMB_e, SSumB, SSumPhi] =...
            Acco(Gmic_compliance_e,  Gmic_compliance_evp,...
            num, Mac_compliance_evp, Interaction_e, Interaction_evp, weight, ...
            dMiceps_0, Micsig, Macsig, dE0, Timeinterval);
        oSumMB_evp(:, :, :, :, inum)  = SSumMB_evp;
        oSumMPhie(:, :, inum)  = SSumMPhie;
        oSumMB_e(:, :, :, :, inum)  = SSumMB_e;
        oSumB(:, :, :, :, inum)  = SSumB;
        oSumPhi(:, :, inum)  = SSumPhi;
    end
    SumMB_evp = sum(oSumMB_evp, 5);
    SumMPhie =  sum(oSumMPhie, 3);
    SumMB_e =  sum(oSumMB_e, 5);
    SumB =  sum(oSumB, 5);
    SumPhi =  sum(oSumPhi, 3);
    % for void
    tempA = Js + doubledotff(inversegeneral(Mac_compliance_evp), Interaction_evp);
    tempAMevp = Mac_compliance_evp + Interaction_evp;
    tempAMe = Mac_compliance_e + Interaction_e;
    temp1 = SumMB_evp + porosity0 * tempAMevp;
    temp2 = SumMB_e + f_p * tempAMe;
    tempO = -doubledotft(Interaction_evp, doubledotft(inversegeneral(Mac_compliance_evp), dE0))...
        - 1/Timeinterval * doubledotft(Interaction_evp, doubledotft(doubledotff(Interaction_e, inversegeneral(Interaction_evp)), Macsig));
    tempAEO = doubledotft(tempA, dE0) + tempO;
    Update_compliance_evp = doubledotff(temp1,  ...
        inversegeneral(SumB));
    Update_dE0 = SumMPhie + porosity0 * tempAEO - doubledotft(Update_compliance_evp, SumPhi);
    Update_compliance_e = doubledotff(temp2,  ...
        inversegeneral(SumB));
    Error1 = abs(norm(tensor2matrix(Update_compliance_evp) - tensor2matrix(Mac_compliance_evp)))/abs(norm(tensor2matrix(Mac_compliance_evp)));
    Error2  = abs(norm(tensor2matrix(Update_compliance_e) - tensor2matrix(Mac_compliance_e)))/abs(norm(tensor2matrix(Mac_compliance_e)));
    Error3 = abs(norm(Update_dE0 - dE0)/norm(dE0));
    Error = max([Error1, Error2]);
    B = zeros(3, 3, 3, 3, length(tnum));
    Phi = zeros(3, 3, length(tnum));
    if Error < 5e-2
        Div = 0;
        Mac_compliance_evp  = Symmetry(Update_compliance_evp);
        Mac_compliance_e  = Symmetry(Update_compliance_e);
        dE0  = Update_dE0;
        Interaction_evp =  doubledotff(doubledotff(inversegeneral(Js - S_evp), S_evp), Mac_compliance_evp);
        for inum = 1 : length(tnum)
            num =  tnum(inum);
            B(:, :, :, :, num) =  doubledotff(inversegeneral(Gmic_compliance_evp(:, :, :, :, num)  +...
                Interaction_evp), (Mac_compliance_evp + Interaction_evp));
            de0 = -1/Timeinterval * doubledotft(Gmic_compliance_e(:, :, :, :, num), Micsig(:, :, num)) + dMiceps_0(:, :, num);
            temp1 = dE0 - de0 -1/Timeinterval * doubledotft(Interaction_e, (Macsig - Micsig(:, :, num)));
            Phi(:, :, num)= doubledotft(inversegeneral(Gmic_compliance_evp(:, :, :, :, num)  +...
                Interaction_evp), temp1);
        end
        break
    else
        j = j + 1;
        if j == j_iter
            Div = 1;
        else
            Div = 0;
        end
        Mac_compliance_evp  = Symmetry(Update_compliance_evp);
        Mac_compliance_e  = Symmetry(Update_compliance_e);
        Mac_compliance_vp  = Mac_compliance_evp - Mac_compliance_e/Timeinterval;
        dE0  = Update_dE0;
        Interaction_evp =  doubledotff(doubledotff(inversegeneral(Js - S_evp), S_evp), Mac_compliance_evp);
        for inum = 1 : length(tnum)
            num =  tnum(inum);
            B(:, :, :, :, num) =  doubledotff(inversegeneral(Gmic_compliance_evp(:, :, :, :, num)  +...
                Interaction_evp), (Mac_compliance_evp + Interaction_evp));
            de0 = -1/Timeinterval * doubledotft(Gmic_compliance_e(:, :, :, :, num), Micsig(:, :, num)) + dMiceps_0(:, :, num);
            temp1 = dE0 - de0 -1/Timeinterval * doubledotft(Interaction_e, (Macsig - Micsig(:, :, num)));
            Phi(:, :, num)= doubledotft(inversegeneral(Gmic_compliance_evp(:, :, :, :, num)  +...
                Interaction_evp), temp1);
        end
    end
end
end

