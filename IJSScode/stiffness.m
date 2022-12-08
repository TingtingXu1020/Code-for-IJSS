function stiffness = stiffness(Micdamage, nvector, num, Identity)
global mu lambda 
                   stiffness = zeros(3, 3, 3, 3);
                   OMEGA_ijkl = zeros(3, 3, 3, 3);
                   NU_ijkl= zeros(3, 3, 3, 3);
                   
                   for i = 1:3
                       for j = 1:3
                           for k = 1:3
                               for l = 1:3
                                   OMEGA_ijkl(i,j,k,l,num) = KronD(i,k)*nvector(num,j)*nvector(num,l)...
                                       +KronD(i,l)*nvector(num,j)*nvector(num,k)...
                                       +KronD(j,k)*nvector(num,i)*nvector(num,l)...
                                       +KronD(j,l)*nvector(num,i)*nvector(num,k);
                                   NU_ijkl(i,j,k,l,num) = nvector(num,i)*nvector(num,j)*nvector(num,k)*nvector(num,l);
                               end
                           end
                       end
                   end
                  
                
                   for i = 1:3
                       for j = 1:3
                           for k = 1:3
                               for l = 1:3
                                   stiffness(i, j, k, l) = lambda*KronD(i,j)*KronD(k,l)+2*mu*Identity(i,j,k,l)...
                                       -Micdamage*mu*OMEGA_ijkl(i,j,k,l,num)...
                                       +(Micdamage)^2*(lambda+2*mu)*NU_ijkl(i,j,k,l,num)...
                                       -Micdamage*lambda*(KronD(i,j)*(nvector(num,k)*nvector(num,l))+(nvector(num,i)*nvector(num,j))*KronD(k,l));
                               end
                           end
                       end
                   end
end