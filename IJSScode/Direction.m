a = 2; b = 2; c = 10;
method = 1;
[angle, rho, R] = Generation(a, b, c, method);
cellnumberA = a*b*c;
cellnumberB = cellnumberA;
cellnumber = cellnumberA + cellnumberB;
weightA = weight;
weightB = weight;
Q(:, :, 1:cellnumberA) = R;
Q(:, :, (cellnumberA + 1) : cellnumber) = R;
weight(1:cellnumberA) = rho;
weight((cellnumberA + 1) : cellnumber) = rho;
%% Bavzant method 42 points for vector
% BETA=33.2699078510*pi/180;
% T1 = cos(pi/4)*cos(pi/2-BETA);
% T2 = cos(BETA);
% n_cosines = [1 0 0;
%            0 1 0;
%            0 0 1;
%            sqrt(2)/2   sqrt(2)/2  0;
%            sqrt(2)/2  -sqrt(2)/2  0;
%            sqrt(2)/2   0          sqrt(2)/2;
%            sqrt(2)/2   0         -sqrt(2)/2;
%            0           sqrt(2)/2  sqrt(2)/2;
%            0           sqrt(2)/2 -sqrt(2)/2;
%            T1          T1         T2;
%            T1          T1        -T2;
%            T1         -T1         T2;
%            T1         -T1        -T2;
%            T1          T2         T1;
%            T1          T2        -T1;
%            T1         -T2         T1;
%            T1         -T2        -T1;
%            T2          T1         T1;
%            T2          T1        -T1;
%            T2         -T1         T1;
%            T2         -T1        -T1;];
% n_cosines(22:42,1:3) = -n_cosines(1:21,1:3);
% n_weight(1:3) = 0.0265214244093;
% n_weight(4:9) = 0.0199301476312;
% n_weight(10:21) = 0.0250712367487;
% n_weight(22:42) = n_weight(1:21);
% temp_1 = acos(n_cosines(:, 3));
% if abs(n_cosines(:, 3)) ~= 1
%     temp_sin = sqrt(1 - n_cosines(:, 3).^2);
% else
%     temp_sin = 1;
% end
% temp_2 = acos(n_cosines(:, 1).*temp_sin.^(-1));
% n_psi = ceil(cellnumber/42);
% temp_3 = 0 : pi/(n_psi-1) : pi;
% angletheta = zeros(42, n_psi);
% anglephi = zeros(42, n_psi);
% anglepsi = zeros(42, n_psi);
% weight1 = zeros(42, n_psi);
% weight2 = zeros(42, n_psi);
% for i = 1 : 42
%     for j = 1 : n_psi
%     angletheta(i, j) = temp_1(i);
%     anglephi(i, j) = temp_2(i);
%     anglepsi(i, j) = temp_3(j);
%     weight1(i, j) = n_weight(i);
%     weight2(i, j) = 1/n_psi;
%     end
% end
% angle (:,1) = reshape(angletheta, 1, []);
% angle (:,2) = reshape(anglephi, 1, []);
% angle (:,3) = reshape(anglepsi, 1, []);
% weight = reshape(weight1, 1, []).*reshape(weight2, 1, [])/2;
%end