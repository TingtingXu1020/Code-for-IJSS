function [Jd, Js, Ja, Jm] = Identity()
% FourIdentity.m
% 4th-order identity tensors
%--------------------------------------------------------------------------

%  allocate J, Jm
   J  = zeros(3,3,3,3);
   J1 = zeros(3,3,3,3);
   Jm = zeros(3,3,3,3);

%  Eqn(3) in Jiang, 2014
   for i = 1:3
       for j = 1:3
           for k = 1:3
               for l = 1:3
                   J(i,j,k,l)  = KronD(i, k)*KronD(j, l);
                   J1(i,j,k,l) = KronD(j, k)*KronD(i, l);
                   Jm(i,j,k,l) = KronD(i, j)*KronD(k, l)/3;
               end
           end
       end
   end
   Js = 0.5*(J + J1);
   Ja = 0.5*(J - J1);
   Jd = Js - Jm;   
end