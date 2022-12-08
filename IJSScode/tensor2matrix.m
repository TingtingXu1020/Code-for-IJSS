%  4th order tensor which has minor symmetry, Mandel notation (Ref, Lee, 2019)

function R = tensor2matrix(F)
    F_matrix(1,1) = F(1,1,1,1);
    F_matrix(1,2) = F(1,1,2,2);
    F_matrix(1,3) = F(1,1,3,3);
    F_matrix(1,4) = sqrt(2)*F(1,1,2,3);
    F_matrix(1,5) = sqrt(2)*F(1,1,3,1);
    F_matrix(1,6) = sqrt(2)*F(1,1,1,2);

    F_matrix(2,1) = F(2,2,1,1);
    F_matrix(2,2) = F(2,2,2,2);
    F_matrix(2,3) = F(2,2,3,3);
    F_matrix(2,4) = sqrt(2)*F(2,2,2,3);
    F_matrix(2,5) = sqrt(2)*F(2,2,3,1);
    F_matrix(2,6) = sqrt(2)*F(2,2,1,2);
    
    F_matrix(3,1) = F(3,3,1,1);
    F_matrix(3,2) = F(3,3,2,2);
    F_matrix(3,3) = F(3,3,3,3);
    F_matrix(3,4) = sqrt(2)*F(3,3,2,3);
    F_matrix(3,5) = sqrt(2)*F(3,3,3,1);
    F_matrix(3,6) = sqrt(2)*F(3,3,1,2);
    
    F_matrix(4,1) = sqrt(2)*F(2,3,1,1);
    F_matrix(4,2) = sqrt(2)*F(2,3,2,2);
    F_matrix(4,3) = sqrt(2)*F(2,3,3,3);
    F_matrix(4,4) = 2*F(2,3,2,3);
    F_matrix(4,5) = 2*F(2,3,3,1);
    F_matrix(4,6) = 2*F(2,3,1,2);
    
    F_matrix(5,1) = sqrt(2)*F(3,1,1,1);
    F_matrix(5,2) = sqrt(2)*F(3,1,2,2);
    F_matrix(5,3) = sqrt(2)*F(3,1,3,3);
    F_matrix(5,4) = 2*F(3,1,2,3);
    F_matrix(5,5) = 2*F(3,1,3,1);
    F_matrix(5,6) = 2*F(3,1,1,2);
    
    F_matrix(6,1) = sqrt(2)*F(1,2,1,1);
    F_matrix(6,2) = sqrt(2)*F(1,2,2,2);
    F_matrix(6,3) = sqrt(2)*F(1,2,3,3);
    F_matrix(6,4) = 2*F(1,2,2,3);
    F_matrix(6,5) = 2*F(1,2,3,1);
    F_matrix(6,6) = 2*F(1,2,1,2);
    
  R = F_matrix;
end
            
        
