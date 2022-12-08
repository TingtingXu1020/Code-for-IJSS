function HScompliance = HScompliance(K, mu, phi, nu)
% Finite Eshelby tensor
alpha = 0; % 0 denotes traction (Neumann) B.C., while 1 denotes displacement (Dirichlet) B.C.
yu_f =  phi * (1 - phi^(2/3))/10/(1-nu)/(7-10*nu);
yt_f =  4 * phi * (1 - phi^(2/3))/10/(1-nu)/(7+5*nu);
ds1d = (1 + nu) / 3/ (1-nu);
ds1n = (1 + nu) / 3/ (1-nu);
ds2inf = 2*(4-5*nu) / 15 / (1-nu);
ds2d =  ds2inf  - 21*yu_f * (1 - phi^(2/3))/(1-phi);
ds2n =  ds2inf + 21 * yt_f *(1- phi^(2/3))/(1-phi);
ds1f =  alpha*ds1d + (1-alpha)*ds1n;
ds2f =  alpha*ds2d + (1-alpha)*ds2n;
K_hom =   K - phi*K*(1-(1 - phi)*ds1f)^(-1);
mu_hom =   mu - phi*mu*(1-(1- phi)*ds2f)^(-1);

% Compliance tensor
E_hom = 9*K_hom*mu_hom/(3*K_hom+mu_hom);
nu_hom = (3*K_hom - 2*mu_hom)/(2*(3*K_hom + mu_hom));
HScompliance = zeros(3, 3, 3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l=1 : 3
                    HScompliance(i, j, k, l) = HScompliance(i, j, k, l)...
                        + (-1)*nu_hom /E_hom *KronD(i,j)*KronD(k,l) + (1 + nu_hom)/2/E_hom*(KronD(i,k)*KronD(j,l)+KronD(i,l)*KronD(j,k));
                end
            end
        end
    end
end