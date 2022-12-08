nu = 0.1;
E= 1E5;
mu = E/(2+2*nu);
Mac_compliance(1:3,1:3,1:3,1:3) = 0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l=1:3
                Mac_compliance(i,j,k,l) = Mac_compliance(i,j,k,l) +...
                    (-1)*nu/E*KronD(i,j)*KronD(k,l) + (1+nu)/2/E*(KronD(i,k)*KronD(j,l)+KronD(i,l)*KronD(j,k));
            end
        end
    end
end

Stiffness = inversegeneral(Mac_compliance);

a = [1;1;1e-15];
N = 16;
M = 4;
Hill = Hill_P(a, N, M, Stiffness);
N_ELE = N*M;
[Alp,Bet,w_pq] = GaussGGLQ(N,M);
Zeta(1:3,1:N_ELE) = 0;
omega(1:N_ELE) = 0;
TGreen(1:3,1:3,1:3,1:3,1:N_ELE) = 0;
Result(1:3,1:3,1:3,1:3,1:N_ELE) = 0;
Hill(1:3,1:3,1:3,1:3) = 0;
Xi(1:3,1:N_ELE) = 0;
K (1:3,1:3,1:N_ELE) = 0;
N_bar (1:3,1:3,1:N_ELE) = 0;
D (1:N_ELE) = 0;

for n_ele = 1:N_ELE
    omega(n_ele) = Alp(n_ele);
    Zeta(3,n_ele) = Bet(n_ele);
    Zeta(1,n_ele) =(1-Zeta(3,n_ele)^2)^(1/2)*cos(omega(n_ele));
    Zeta(2,n_ele) =(1-Zeta(3,n_ele)^2)^(1/2)*sin(omega(n_ele));
    
    for i = 1:3
        Xi(i,n_ele) = Zeta(i,n_ele)/a(i);
    end
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    K (i,k,n_ele) = K (i,k,n_ele) + Stiffness (i,j,k,l) * Xi(j,n_ele) * Xi(l,n_ele);
                end
            end
        end
    end
    
    for m = 1:3
        for n = 1:3
            for l = 1:3
                D(n_ele) =  D(n_ele) + LeviCivita([m,n,l]) * K(m,1,n_ele) * K(n,2,n_ele) * K(l,3,n_ele);
            end
        end
    end
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        for n =1:3
                            N_bar(i,j,n_ele) = N_bar(i,j,n_ele) + 1/2 * LeviCivita([i,k,l]) * LeviCivita([j,m,n]) * K(k,m,n_ele) * K(l,n,n_ele);
                        end
                    end
                end
            end
        end
    end
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    TGreen(i,j,k,l,n_ele) = TGreen(i,j,k,l,n_ele)  + Xi(k,n_ele) * Xi(l,n_ele) *N_bar(i,j,n_ele)/D(n_ele);
                end
            end
        end
    end
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3                  
                    Result(i,j,k,l,n_ele) = Result(i,j,k,l,n_ele) + (TGreen(i,l,j,k,n_ele) + TGreen(j,l,i,k,n_ele));    
                end
            end
        end
    end
 
    Hill = Hill + 1/(8*pi) * (1/2) * (2*pi-0) * (Result(:,:,:,:,n_ele)*w_pq(n_ele));
end
T = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3                  
                    T(i, j, k, l) = T(i,j,k,l) + 1/4 * (Hill(i, j, k, l) + Hill(i, j, l, k) + Hill(j, i, k, l) + Hill(j, i, l, k));    
                end
            end
        end
    end
M_1 = doubledotff(T, Stiffness);

tensor2matrix(M_1)

S(1:3,1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l=1:3
                S(i,j,k,l) = (5*nu-1)/(15*(1-nu))*KronD(i,j)*KronD(k,l) + (4-5*nu)/(15*(1-nu))*(KronD(i,k)*KronD(j,l)+KronD(i,l)*KronD(j,k));
            end
        end
    end
end
[Jd, Js, Ja, Jm] = Identity();
tensor2matrix(S)

Interaction = inversegeneral(doubledotff(doubledotff(inversegeneral(Js - S),S), Mac_compliance));


tensor2matrix(Interaction);
Ltensor_intact(1:3,1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l=1:3
                Ltensor_intact(i,j,k,l) = mu/(4-5*nu)*((3-5*nu)*KronD(i,j)*KronD(k,l)+(7-5*nu)/2*(KronD(i,k)*KronD(j,l)+KronD(i,l)*KronD(j,k)));
            end
        end
    end
end

tensor2matrix(Ltensor_intact);