function [NdMicsig] ...
    = local2REV(num, dMaceps, Interaction_vp, Macsig, Interaction_e, ...
    dMacsig, Gmic_compliance_vp, Micsig, Gmic_compliance_e, sigma_p, Timeinterval)
% Use the stress from last step 
tem_b = dMaceps - doubledotft(Interaction_vp, (Micsig(:,:,num) - Macsig)) + doubledotft(Interaction_e, dMacsig)...
      - doubledotft( Gmic_compliance_vp(:,:,:,:,num) , Micsig(:,:,num) - sigma_p (:,:,num));
tem_A = Interaction_e  + Gmic_compliance_e(:,:,:,:,num) + (Interaction_vp  + Gmic_compliance_vp(:,:,:,:,num))*Timeinterval;
% Use the stress this step
NdMicsig = doubledotft(inversegeneral(tem_A), tem_b);
end