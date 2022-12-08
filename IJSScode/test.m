inc = 5;
Titer = 1;
j_vp = 1;
 while j_vp < 100
                Accomodation_vp(1:3,1:3,1:3,1:3,1:cellnumber,inc,j_vp,Titer)=0;
                Interaction_vp(1:3,1:3,1:3,1:3,inc,j_vp,Titer)=0;
                Sum_vp(1:3,1:3,1:3,1:3) =0;
                Sum_acco_vp(1:3,1:3,1:3,1:3) =0;
%                 if det(tensor2matrix(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer))) < 10^(-30)
%                     Mod_compliance_vp = tensor2matrix(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer)) + max(max(max(max(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer))))) *eye(6)*10^(-2);
%                     Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer) = matrix2tensor(Mod_compliance_vp);
%                 end
                S_E_vp(:,:,:,:,inc,j_vp,Titer) = Eshelby(Cell_R*[1;1;1], 16,16, inversegeneral(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer)));
                % Temp(:,:,:,:,inc,j_vp,Titer) = doubledotff(doubledotff(inversegeneral(Identity - S_E_vp(:,:,:,:,inc,j_vp,Titer)),S_E_vp(:,:,:,:,inc,j_vp,Titer)), Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer));
                Temp_vp(:,:,:,:,inc,j_vp,Titer) = doubledotff(inversegeneral(inversegeneral(S_E_vp(:,:,:,:,inc,j_vp,Titer)) - Identity) ,  Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer));
                Interaction_vp(:,:,:,:,inc,j_vp,Titer) = Temp_vp(:,:,:,:,inc,j_vp,Titer);
                Ltensor_vp(:,:,:,:,inc,j_vp,Titer) = inversegeneral(Interaction_vp(:,:,:,:,inc,j_vp,Titer));
                for num = 1:cellnumber
                    Accomodation_vp(:,:,:,:,num,inc,j_vp,Titer) = doubledotff(inversegeneral(Gmic_compliance_vp(:,:,:,:,num,inc)  +...
                        Interaction_vp(:,:,:,:,inc,j_vp,Titer)),(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer) + Interaction_vp(:,:,:,:,inc,j_vp,Titer)));
                    Sum_vp  = Sum_vp + doubledotff(Gmic_compliance_vp(:,:,:,:,num,inc),(Accomodation_vp(:,:,:,:,num,inc,j_vp,Titer)));
                    Sum_acco_vp = Sum_acco_vp + Accomodation_vp(:,:,:,:,num,inc,j_vp,Titer);
                end
                
                if abs(norm(tensor2matrix(doubledotff(Sum_vp/cellnumber,  inversegeneral(Sum_acco_vp/cellnumber))) - tensor2matrix(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer)))) <= 10^(-2)*abs(norm(tensor2matrix(Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer))))
                    j_vp = j_vp + 1;
                    Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer)  = doubledotff(Sum_vp/cellnumber,  inversegeneral(Sum_acco_vp/cellnumber));
                    break
                else
                    j_vp = j_vp + 1;
                    Mac_compliance_vp(:,:,:,:,inc,j_vp,Titer)  = doubledotff(Sum_vp/cellnumber,  inversegeneral(Sum_acco_vp/cellnumber));
                end
            end
