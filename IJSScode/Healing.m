function [dvoid_r, Lmic_compliance_c]  = Healing(LMicsig, void_r, Cell_R, scale_up, R, T, P_l)
% Parameters
% D = 0.0013;   % mm^2/s
% S = 2e-6;  % mm
Omega = 2.7e4; % mm^3/mol
% Cpore = 6.48e-6;  % mol/mm^3
Z_star = 8.16e-13; %mm^3/s
% Initial set
sigma_n = diag((LMicsig + void_r.^3/Cell_R^3 * P_l *eye(3))/(1 - void_r.^3/Cell_R^3));
Gps = -(Z_star * Omega^2*(Cell_R^2 - void_r.^2))./...
    (R*T*((3*Cell_R^4)/4 - Cell_R^2*void_r.^2 + void_r.^4/4 - Cell_R^4*log(Cell_R) + Cell_R^4*log(void_r))); % convergence rate factor
Ac = 2 * pi * (Cell_R^2 - void_r.^2); 
V = -sigma_n.* Gps; % compressive stress
V(V<0) = 0;
Vv = Ac.*V; % dissolved volume rate
Vt = sum(Vv);
if  void_r < 1e-15
    dvoid_r = 0;
else
    Isurface = 4 * pi * void_r .^2;
    dvoid_r = -Vt/Isurface;
end
if max (dvoid_r) > 0
    error('error')
end

healing_matrix = scale_up*1/Cell_R*(diag(Gps) * diag(sign(V)));
Lmic_compliance_c = zeros(3,3,3,3);

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l=1:3
                Lmic_compliance_c(i, j, k, l)= Lmic_compliance_c(i, j, k, l) + ...
                    healing_matrix(i, j) * (KronD(i, k) * KronD(j,l) * KronD(i,j) * KronD(k,l))  ;
            end
        end
    end
end

end