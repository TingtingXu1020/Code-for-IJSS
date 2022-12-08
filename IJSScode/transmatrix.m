function [Q] = transmatrix(theta1, phi1, theta2, phi2)
%function [Q] = transmatrix(theta, phi, psi)
% R1 = [cos(psi) sin(psi) 0
%     -sin(psi) cos(psi) 0
%     0 0 1];
% 
% R2 = [cos(theta) 0 sin(theta)
%     0 1 0
%     -sin(theta) 0 cos(theta)];
% 
% R3 = [cos(phi) sin(phi) 0
%     -sin(phi) cos(phi) 0
%     0 0 1];
% 
% Q = R1*R2*R3;

e1 = [cos(theta1)*sin(phi1); sin(theta1)*sin(phi1); cos(phi1)];
if phi1 == pi/2
     e2 = [-sin(theta1)*sin(phi2); cos(theta1)*sin(phi2); cos(phi2)];
else
   e2 = [cos(theta2)*cos(atan(cos(theta1-theta2)*tan(phi1))); ...
       sin(theta2)*cos(atan(cos(theta1-theta2)*tan(phi1))); -sin(atan(cos(theta1-theta2)*tan(phi1)))];
end
e3 = cross(e1, e2);
Q = [transpose(e1); transpose(e2); transpose(e3)];

end