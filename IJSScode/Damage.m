function [Ldamage_matrix, Lmic_compliance_e , NLmic_compliance_e, LMicsig, Nlmicdamage] ...
    = Damage( dMiceps, Micsig,  lambda, mu,  Q, lmicdamage, sigtensile, eta, Timeinterval, d_c)
% Projection on contact planes and Damage Criterion Check
dleps =  Q * dMiceps * transpose(Q);
lsig = Q * Micsig * transpose(Q);
Tr_lsig =  lsig + doubledotft(Damage_stiffness(lmicdamage, lambda, mu), dleps)*Timeinterval;
Maxiter = 1E3;
% Damage elastic trial guess (i.e., no damage evolution)
%% consider three boundary interfaces
sigma_n = diag(Tr_lsig);
fd = diag(Tr_lsig) -  sigtensile*ones(3,1);
%% Solve the consistent condition
% Calculate the dimension of the damage active system
iter = 0;
z = find(fd > 0 & sigma_n>0);
Nlmicdamage = lmicdamage;
if max (fd) <  1E-6 || isempty(z)% No damage
    dlmicdamage = zeros(3,1);
    Nlmicdamage = lmicdamage + dlmicdamage;
else % Consistent condition for damage increment
    while max(fd) > 1E-6 && iter < Maxiter
        uplmicdamage = Activedamage(z, lmicdamage, lambda, mu, dleps, fd, lsig,  sigtensile, eta, d_c);
        % Recalculate damage function regarding the maximum damage systems
        check_lsig =  lsig + doubledotft(Damage_stiffness(uplmicdamage, lambda, mu), dleps);
        xi = uplmicdamage/d_c;
        dhardening = (1 + eta * xi./(1 + xi.^2));
        fd = diag(check_lsig) -   sigtensile*dhardening;
        z = unique(sort([z; find(fd > 1E-6)]));
        iter = iter + 1;
        Nlmicdamage = uplmicdamage;
        if abs(Nlmicdamage(Nlmicdamage ~= 0) - 1) < 1E-5
            break
        end
        %         if iter == Maxiter
        %             disp('Damage part is diverse');
        %         end
    end
end
% output
Nlmicdamage(Nlmicdamage > 1) = 1 - 1E-5;
Ldamage_matrix =  [Nlmicdamage(1), 0, 0; 0, Nlmicdamage(2), 0; 0, 0, Nlmicdamage(3)];
Lmic_compliance_e =  inversegeneral(Damage_stiffness(lmicdamage, lambda, mu));
NLmic_compliance_e = inversegeneral(Damage_stiffness(Nlmicdamage, lambda, mu));
LMicsig = lsig + doubledotft(Damage_stiffness(Nlmicdamage, lambda, mu), dleps);
end