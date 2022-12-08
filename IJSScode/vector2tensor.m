% Applied for symmetric 2nd order tensor, for example, stress and strain
% tensor (from 3-by-3 matrix to 6-by-1 vector) Mandel notation

function tensor = vector2tensor(v)
tensor(1, 1) = v(1);
tensor(2, 2) = v(2);
tensor(3, 3) = v(3);
tensor(2, 3) = v(4)/sqrt(2);
tensor(3, 1) = v(5)/sqrt(2);
tensor(1, 2) = v(6)/sqrt(2);

tensor(3, 2) = tensor(2, 3);
tensor(1, 3) = tensor(3, 1);
tensor(2, 1) = tensor(1, 2);
end