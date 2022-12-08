function uplmicdamage = Activedamage(z, lmicdamage, lambda, mu, dleps, fd, lsig,  sigtensile, eta, d_c)
iter = 0;
Maxiter = 1E3;
uplmicdamage = lmicdamage;
while iter <   Maxiter% Local Newton iteration
    Jd = Jacob(lmicdamage, lambda, mu, dleps, sigtensile, eta, d_c);
    dlmicdamage = zeros(3, 1);
    if length(z) == 1
        dlmicdamage(z) =  -fd(z) * Jd(z,z)^(-1);
    elseif length(z) == 2
        Jd_2D = [Jd(z(1), z(1)), Jd(z(1), z(2))
            Jd(z(2), z(1)), Jd(z(2), z(2))];
        if rcond(Jd_2D) < 1E-15
            ddamage_2D = -pinv(Jd_2D) * [fd(z(1)); fd(z(2))];
        else
            ddamage_2D = -Jd_2D\ [fd(z(1)); fd(z(2))];
        end
        dlmicdamage(z(1)) = ddamage_2D(1);
        dlmicdamage(z(2)) = ddamage_2D(2);
    elseif length(z) == 3
        if rcond(Jd) < 1E-15
            dlmicdamage =  - (pinv(Jd)  * [fd(1); fd(2); fd(3)]);
        else
            dlmicdamage =  - (Jd\[fd(1); fd(2); fd(3)]);
        end
    end
    % Drop the inactive damage system
    if isempty(find(dlmicdamage > 0, 1))
        dlmicdamage = zeros(3,1);
        uplmicdamage = uplmicdamage + dlmicdamage;
        break
    end
    while min(dlmicdamage) < 0
        if isempty(find(dlmicdamage > 0, 1))
            dlmicdamage = zeros(3,1);
            break
        end
        z = setdiff(z, find(dlmicdamage < 0));
        dlmicdamage(dlmicdamage < 0) = 0;
    end
    % Solve the damage function for current active system    
    uplmicdamage = uplmicdamage + dlmicdamage;
    if abs(uplmicdamage(uplmicdamage ~= 0) - 1) < 1E-5
        break
    end
    uplmicdamage(uplmicdamage > 1) = 1 - 1E-5;
    check_lsig =  lsig + doubledotft(Damage_stiffness(uplmicdamage, lambda, mu), dleps);
    xi = uplmicdamage/d_c;
    dhardening = (1 + eta * xi./(1 + xi.^2));
    fd = diag(check_lsig) -   sigtensile*dhardening;

    %  Termination condition
    if max(fd(z)) < 1E-6
        break
    end
    iter = iter + 1;
%     if iter == Maxiter
%         disp('Local Damage is diverse');
%     end
end
end