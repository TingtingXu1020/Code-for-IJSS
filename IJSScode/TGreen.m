function T      = TGreen(R, C_array, Alp, Bet, ww)
% TGreen.m
% 4th-order Green tensor T
% 
% Input:  x, three semi-axes of the inclusion, 3*1 matrix;
%         C, stiffness of the matrix, 1*81 matrix;
%         Alp#, Bet#, ww#, nodes and weights;
%       
% Output: T, 4th-order tensor, 3*3*3*3 matrix;
%--------------------------------------------------------------------------
   T             = zeros(3,3,3,3);
   for i = 1:3
       for j = 1:3
           for k = i:3
               for l = j:3
                       T(i,j,k,l) = GGLQ([i,j,k,l], R, 0, 2*pi, 0, pi, Alp, Bet, ww, C_array);
                     %  T(i,j,k,l) = GLeQ([i,j,k,l], R, Alp, Bet, ww, C_array);
                       T(k,j,i,l) = T(i,j,k,l);
                       T(i,l,k,j) = T(i,j,k,l);
                       T(k,l,i,j) = T(i,j,k,l);
               end
           end
       end
  end   
end
