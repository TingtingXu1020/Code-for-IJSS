function R = outerproducttt(A,B)
R = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l =1:3
                R(i, j, k, l) = R(i, j, k, l) + A(i, j)*B(k, l);
            end
        end
    end
end
end
