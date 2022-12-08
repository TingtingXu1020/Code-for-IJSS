function [angle, weight, R, number] = Generation(a, b, c, method)
switch method
    case 1
        %% Jiang's method
        if a == b && b ==c % Randomly linear combination
            number = a;
            theta1= 0 : 2*pi/(a-1) : 2*pi;
            theta1 = theta1(randperm(a));
            phi1= acos( -1 : 2/(b-1) : 1);
            phi1 = phi1(randperm(b));
            psi= 0 : pi/(c-1) : pi;
            psi = psi(randperm(c));
            phi2= acos( -1 : 2/(c-1) : 1);
            phi2 = phi2(randperm(c));
            theta2 = theta1 + pi/2 + atan(cos(phi1) .* tan(psi));
            
            angle (:,1) = theta1;
            angle (:,2) = phi1;
            angle (:,3) = theta2;
            angle (:, 4) = phi2;
            weight(1 : number) = 1/number;
        else % Cubic combination
            
            number = a*b*c;
            if a == 1
                theta1 = pi/4;
            else
            theta1= 0 : 2*pi/(a-1) : 2*pi;
            theta1 = theta1(randperm(a));
            end
            phi1= acos( -1 : 2/(b-1) : 1);
            phi1 = phi1(randperm(b));
            psi= 0 : pi/(c-1) : pi;
            psi = psi(randperm(c));
            phi2= acos( -1 : 2/(c-1) : 1);
            phi2 = phi2(randperm(c));
            
            theta1 = repmat(theta1', b*c, 1);
            phi1 = reshape(repmat(phi1, a, c), a*b*c, 1);
            psi = reshape(repmat(psi, a*b, 1), a*b*c, 1);
            phi2 = reshape(repmat(phi2, a*b, 1), a*b*c, 1);
            
            theta2 = theta1 + pi/2 + atan(cos(phi1) .* tan(psi));
            
            angle (:,1) = theta1;
            angle (:,2) = phi1;
            angle (:,3) = theta2;
            angle (:, 4) = phi2;
            weight(1 : number) = 1/number;
        end
        
        % Transformation matrix
        R(1:3, 1:3, 1:number) = 0;
        for num = 1:number
            e1 = [cos(theta1(num))*sin(phi1(num)); sin(theta1(num))*sin(phi1(num)); cos(phi1(num))];
            if phi1 == pi/2
                e2 = [-sin(theta1(num))*sin(phi2(num)); cos(theta1(num))*sin(phi2(num)); cos(phi2(num))];
            else
                e2 = [cos(theta2(num))*cos(atan(cos(theta1(num)-theta2(num))*tan(phi1(num)))); ...
                    sin(theta2(num))*cos(atan(cos(theta1(num)-theta2(num))*tan(phi1(num)))); -sin(atan(cos(theta1(num)-theta2(num))*tan(phi1(num))))];
            end
            e3 = cross(e1, e2);
            R(:, :, num) = [transpose(e1); transpose(e2); transpose(e3)];
        end
    case 2
        %% Pouya's method
        number = a*b*c;
        if a == 1
            psi = pi/4;
        else
        psi = 0 : 2*pi/(a-1) : 2*pi;
        end
        psi = psi(randperm(a));
        theta= acos( -1 : 2/(b-1) : 1);
        theta = theta(randperm(b));
        phi= 0 : 2*pi/(c-1) : 2*pi;
        phi = phi(randperm(c));
        
        psi = repmat(psi', b*c, 1);
        theta = reshape(repmat(theta, a, c), a*b*c, 1);
        phi = reshape(repmat(phi, a*b, 1), a*b*c, 1);
        
        
        angle (:,1) = psi ;
        angle (:,2) = theta;
        angle (:,3) = phi;
        weight(1 : number) = 1/number;
        
        % Transformation matrix
        R(1:3, 1:3, 1:number) = 0;
        for num = 1 : number
            R1 = [cos(psi(num)) sin(psi(num)) 0
                -sin(psi(num)) cos(psi(num)) 0
                0 0 1];
            
            R2 = [cos(theta(num)) 0 sin(theta(num))
                0 1 0
                -sin(theta(num)) 0 cos(theta(num))];
            
            R3 = [cos(phi(num)) -sin(phi(num)) 0
                sin(phi(num)) cos(phi(num)) 0
                0 0 1];
            
            
            R(:, :, num) = R1*R2*R3;
        end
end
end