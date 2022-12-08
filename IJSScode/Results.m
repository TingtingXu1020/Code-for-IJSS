load('DoS1Maceps_save.mat')
load('DoS1Maceps_v.mat')
load('DoS1Macporosity.mat')
load('DoS1Macsig_save.mat')
load('DoS1Micsig.mat')
load('DoS1T.mat')

DOS1Maceps = Maceps_save;
DOS1Maceps_v = Maceps_v;
DOS1T = T_save;


load('DoS0.5Maceps_save.mat')
load('DoS0.5Maceps_v.mat')
load('DoS0.5Macporosity.mat')
load('DoS0.5Macsig_save.mat')
load('DoS0.5Micsig.mat')
load('DoS0.5T.mat')

DOShalfMaceps_v = Maceps_v;


plot(DOS1T, DOS1Maceps_v(1:end-1))
hold on
%plot(DOShalfMaceps_v)