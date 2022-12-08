load('DoS1Cell_R0.2Active2T.mat')
load('DoS1Cell_R0.2Active2Micsig.mat')
load('DoS1Cell_R0.2Active2Macsig_save.mat')
load('DoS1Cell_R0.2Active2Macporosity.mat')
load('DoS1Cell_R0.2Active2Maceps_v.mat')
load('DoS1Cell_R0.2Active2Maceps_save.mat')
load('DoS1Cell_R0.2Active2Mac_compliance_vpsave.mat')
load('DoS1Cell_R0.2Active2Mac_compliance_esave.mat')


for i =  1 : size(Mac_compliance_esave, 5)
    E1111(i) = Mac_compliance_esave(1, 1, 1, 1, i);
    E2222(i) = Mac_compliance_esave(2, 2, 2, 2, i);
    E3333(i) = Mac_compliance_esave(3, 3, 3, 3, i);
    E2323(i) = Mac_compliance_esave(2, 3, 2, 3, i);
    E1313(i) = Mac_compliance_esave(1, 3, 1, 3, i);
    E1212(i) = Mac_compliance_esave(1, 2, 1, 2, i);
end
close all
plot(E1111(1:end))
%%

for i =  1 : size(Mac_compliance_vpsave, 5)
    V1111(i) = Mac_compliance_vpsave(1, 1, 1, 1, i);
    V2222(i) = Mac_compliance_vpsave(2, 2, 2, 2, i);
    V3333(i) = Mac_compliance_vpsave(3, 3, 3, 3, i);
    V2323(i) = Mac_compliance_vpsave(2, 3, 2, 3, i);
    V1313(i) = Mac_compliance_vpsave(1, 3, 1, 3, i);
    V1212(i) = Mac_compliance_vpsave(1, 2, 1, 2, i);
end

close all
plot(V1111(1:100))