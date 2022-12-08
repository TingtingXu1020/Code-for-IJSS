%A is 4th order tensor

function R = normh(A)

R = 0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                R = R + A(i,j,k,l)*A(i,j,k,l);
            end
        end
    end
end

R = sqrt(R);