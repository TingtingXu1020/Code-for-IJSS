function Jd = Jacob(lmicdamage, lambda, mu, dleps, sigtensile, eta, d_c)
lmicdamage(lmicdamage == 1) = 1 - 1E-5;
%% Softening
Jd = [ - dleps(1, 1)*(lambda + 2*mu) - sigtensile*(eta/(d_c*(lmicdamage(1)^2/d_c^2 + 1)) - (2*eta*lmicdamage(1)^2)/(d_c^3*(lmicdamage(1)^2/d_c^2 + 1)^2)) - (dleps(2, 2)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(1))^(1/2)) - (dleps(3, 3)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(1))^(1/2)),                                                                                                                                                                                                                  -(dleps(2, 2)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(2))^(1/2)),                                                                                                                                                                                                                  -(dleps(3, 3)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(3))^(1/2))
                                                                                                                                                                                                                 -(dleps(1, 1)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(1))^(1/2)), - dleps(2, 2)*(lambda + 2*mu) - sigtensile*(eta/(d_c*(lmicdamage(2)^2/d_c^2 + 1)) - (2*eta*lmicdamage(2)^2)/(d_c^3*(lmicdamage(2)^2/d_c^2 + 1)^2)) - (dleps(1, 1)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(2))^(1/2)) - (dleps(3, 3)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(2))^(1/2)),                                                                                                                                                                                                                  -(dleps(3, 3)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(3))^(1/2))
                                                                                                                                                                                                                  -(dleps(1, 1)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(1))^(1/2)),                                                                                                                                                                                                                  -(dleps(2, 2)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(2))^(1/2)), - dleps(3, 3)*(lambda + 2*mu) - sigtensile*(eta/(d_c*(lmicdamage(3)^2/d_c^2 + 1)) - (2*eta*lmicdamage(3)^2)/(d_c^3*(lmicdamage(3)^2/d_c^2 + 1)^2)) - (dleps(1, 1)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(3))^(1/2)) - (dleps(2, 2)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(3))^(1/2))];
 
%% Hardening 
% Jd11 = lmicdamage(1)/(sigtensile*(lmicdamage(1)/sigtensile + 1/eta)^2) - 1/(lmicdamage(1)/sigtensile + 1/eta) - eps(1,1)*(lambda + 2*mu) - ...
%     (eps(2,2)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(1))^(1/2)) - (eps(3,3)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(1))^(1/2));
% Jd12 = -(eps(2,2)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(2))^(1/2));
% Jd13 = -(eps(3,3)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(3))^(1/2));
% 
% Jd21 =-(eps(1,1)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(1))^(1/2));
% Jd22 = lmicdamage(2)/(sigtensile*(lmicdamage(2)/sigtensile + 1/eta)^2) - 1/(lmicdamage(2)/sigtensile + 1/eta) - eps(2,2)*(lambda + 2*mu) - ...
%     (eps(1,1)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(2))^(1/2)) - (eps(3,3)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(2))^(1/2));
% Jd23 = -(eps(3,3)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(3))^(1/2));
%  
% 
% Jd31 = -(eps(1,1)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(1))^(1/2));
% Jd32 = -(eps(2,2)*lambda*(1 - lmicdamage(3))^(1/2))/(2*(1 - lmicdamage(2))^(1/2));
% Jd33 = lmicdamage(3)/(sigtensile*(lmicdamage(3)/sigtensile + 1/eta)^2) - 1/(lmicdamage(3)/sigtensile + 1/eta) - eps(3,3)*(lambda + 2*mu) - ...
%     (eps(1,1)*lambda*(1 - lmicdamage(1))^(1/2))/(2*(1 - lmicdamage(3))^(1/2)) - (eps(2,2)*lambda*(1 - lmicdamage(2))^(1/2))/(2*(1 - lmicdamage(3))^(1/2));
%  

% 
% Jeps(1,1)= -(lambda + 2*mu)*(1 - lmicdamage(1))^(1/2)*(1 - lmicdamage(1))^(1/2);
% Jeps(1,2) = lambda*(1 - lmicdamage(1))^(1/2)*(1 - lmicdamage(2))^(1/2);
% Jeps(1,3) = lambda*(1 - lmicdamage(1))^(1/2)*(1 - lmicdamage(3))^(1/2);
% Jeps(2,1) = lambda*(1 - lmicdamage(1))^(1/2)*(1 - lmicdamage(2))^(1/2);
% Jeps(2,2) =  - (lambda + 2*mu)*(1 - lmicdamage(2))^(1/2)*(1 - lmicdamage(2))^(1/2);
% Jeps(2,3) = lambda*(1 - lmicdamage(2))^(1/2)*(1 - lmicdamage(3))^(1/2);
% Jeps(3,1) = lambda*(1 - lmicdamage(1))^(1/2)*(1 - lmicdamage(3))^(1/2);
% Jeps(3,2) = lambda*(1 - lmicdamage(2))^(1/2)*(1 - lmicdamage(3))^(1/2);
% Jeps(3,3)=   - (lambda + 2*mu)*(1 - lmicdamage(3))^(1/2)*(1 - lmicdamage(3))^(1/2);

end





