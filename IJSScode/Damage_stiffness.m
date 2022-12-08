function C = Damage_stiffness(lmicdamage, lambda, mu)
D = diag(sqrt(ones(3, 1) - lmicdamage));
C = zeros(3, 3, 3, 3);
for i = 1:3
    for j = 1:3
         for k = 1:3
             for l = 1:3
                 C(i, j, k, l) = C(i, j, k, l) + lambda*D(i, j)*D(k, l) + mu*(D(i, k)*D(j, l) + D(i, l) *D(j, k));
             end
         end
    end
end
end