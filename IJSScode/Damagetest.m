% Elastic parameter (Ref. Liang 2006)
E = 1E4; % Young's modulus
nu = 0.31;  % Possion's ratio
K = E/3/(1 - 2*nu);
lambda = E*nu /(1 + nu)/(1 - 2*nu);
mu = E/2/(1+nu);
Timeinterval = 1;
lmicdamage = zeros(3, 1);
% Initial condition
Q = [ -0.6103   -0.7918    0.0234
    0.7862   -0.6091   -0.1048
    0.0973   -0.0456    0.9942];
Micsig = zeros(3, 3);
LMicsig = zeros(3, 3);
dMiceps = [1, 0, 0; 0, 0, 0; 0, 0, 0]*1E-5;
Miceps = zeros(3, 3);
% Intact elastic stiffness
Lmic_stiffness_e = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l=1:3
                Lmic_stiffness_e(i,j,k,l) = Lmic_stiffness_e(i,j,k,l)...
                    + lambda *KronD(i,j)*KronD(k,l) + mu*(KronD(i,k)*KronD(j,l) + KronD(i,l)*KronD(j,k));
            end
        end
    end
end

T = 0;

while norm(Miceps) < 1E-3
    % Elastic trial guess for damage model
    Miceps = Miceps + dMiceps* Timeinterval;
    % Projection on contact planes and Damage Criterion Check
    sigtensile = 3;                 % Tensile strength (MPa)
    eta = 0.5;                      % Larger eta gives smaller numerical difficulty
    ncosine = eye(3);
    lsig = Q * Micsig * transpose(Q);
    leps =  Q * Miceps * transpose(Q);
    dleps = Q * dMiceps * transpose(Q);
    Tr_lsig =  lsig + doubledotft(Lmic_stiffness_e, dleps)*Timeinterval;
    sigma_n(1:3) = 0;
    eps_n(1:3) = 0;
    fd (1:3) = 0;
    for i_point = 1:3
        sigma_n(i_point) = ncosine(i_point, :) * Tr_lsig * transpose(ncosine(i_point, :));
        eps_n(i_point) = ncosine(i_point, :) * leps * transpose(ncosine(i_point, :));
        fd(i_point) = sigma_n(i_point) - ( (sigtensile -  E *eps_n(i_point)  + eta*sigtensile) /eta);
    end
    
    
    % Calculate the dimension of the linear system
    z = find(fd > 0);
    if max (fd) < 10^(-6) % No damage
        lmicdamage = zeros(3, 1);
    else % Consistent condition for damage increment
        iter = 1;
         while iter < 1E3
            [Jd, Jeps] = Jacob(lmicdamage, lambda, mu, leps, E, eta);
            feps = Jeps*[dleps(1, 1); dleps(2, 2); dleps(3, 3)];
            dlmicdamage = zeros(3,1);
            
            if length(z) == 1
                dlmicdamage(z) =  -feps(z) * Jd(z,z)^(-1);
            elseif length(z) == 2
                Jd_2D = [Jd(z(1), z(1)), Jd(z(1), z(2))
                    Jd(z(2), z(1)), Jd(z(2), z(2))];
                if cond(Jd_2D) > 10^(15)
                    ddamage_2D = -pinv(Jd_2D) * [feps(z(1)); feps(z(2))];
                else
                    ddamage_2D = -Jd_2D\ [feps(z(1)); feps(z(2))];
                end
                dlmicdamage(z(1)) = ddamage_2D(1);
                dlmicdamage(z(2)) = ddamage_2D(2);
            elseif length(z) == 3
                if cond(Jd) > 10^(15)
                    dlmicdamage =  - (pinv(Jd)  * [feps(1); feps(2); feps(3)]);
                else
                    dlmicdamage =  - (Jd\[feps(1); feps(2); feps(3)]);
                end
            end
            
            for i = 1:3
                if dlmicdamage(i) < 0
                    dlmicdamage(i) = 0;
                end
            end
            
            lmicdamage = lmicdamage + dlmicdamage;
            check_lsig =  LMicsig + doubledotft(Damage_stiffness(lmicdamage, lambda, mu), dleps) ;
            for i = 1 : length(z)
                sigma_n(z(i)) = ncosine(z(i), :)  * check_lsig *  transpose(ncosine(z(i), :));
                fd(z(i)) = sigma_n(z(i)) - ( (sigtensile -  E *eps_n(z(i))  + eta*sigtensile) /eta);
            end
            iter = iter + 1;
               
                        % Termination condition
            if max(fd) < 10^(-2)
                break;
            end
        end
    end
    
    
    
    % output
    Ldamage_matrix =  [lmicdamage(1), 0, 0; 0, lmicdamage(2), 0; 0, 0, lmicdamage(3)];
    Lmic_compliance_e = inversegeneral(Damage_stiffness(lmicdamage, lambda, mu));
    LMicsig =  LMicsig + doubledotft(Damage_stiffness(lmicdamage, lambda, mu), dleps) ;
    Micsig = transpose(Q) * LMicsig * Q;
    
    T = T + 1;
    X(T) = Miceps(1, 1);
    y1(T) = Micsig(1, 1);
    y2(T) = Micsig(2, 2);
    y3(T) = Micsig(3, 3);
    D1(T) = lmicdamage(1);
    D2(T) = lmicdamage(2);
    D3(T) = lmicdamage(3);
end

plot(X, y1);
plot(X, D1);