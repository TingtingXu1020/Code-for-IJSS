% Applied for symmetric 2nd order tensor, for example, stress and strain
% tensor (from 3-by-3 matrix to 6-by-1 vector) Mandel notation

function vector = tensor2vector(M)
vector = zeros(6, 1);
vector(1) = M(1, 1);
vector(2) = M(2, 2);
vector(3) = M(3, 3);
vector(4) = M(2, 3);
vector(5) = M(3, 1);
vector(6) = M(1, 2);
end
