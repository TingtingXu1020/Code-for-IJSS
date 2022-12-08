%Both A and B are 2nd order tensors

function R = otimestt(A,B)

R(1:3,1:3,1:3,1:3) = 0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                R(i,j,k,l) = A(i,j)*B(k,l);
            end
        end
    end
end




