function S = Symmetry(C)
S(1:3, 1:3, 1:3, 1:3) = 0;
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l=1:3
                          S(i,j,k,l) = S(i,j,k,l) + (C(i,j,k,l) + C(k,l,i,j))/2; 
                        end
                    end
                end
            end   
end