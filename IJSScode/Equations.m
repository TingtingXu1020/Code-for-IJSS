syms A B r Omega Vps R T D C S sigma
Ac = pi * (A^2 - r^2);
Ap = 2 * pi * r * S;
ddotdeltar = 4 * Vps^2 *Ac^2/Omega^2/Ap * R * T/D/C;
%simplify(ddotdeltar)
%dotdeltar = simplify(int(ddotdeltar, r));
dotdeltar = (R*T*Vps^2*pi*(A^4 - 4*A^2*A^2 + 4*A^4*log(A)))/(2*C*D*Omega^2*S)...
    - (R*T*Vps^2*pi*(B^4 - 4*A^2*B^2 + 4*A^4*log(B)))/(2*C*D*Omega^2*S);
Vpssolution = solve(simplify(dotdeltar) == 2* pi * (A^2 - B^2)*Vps*sigma, Vps );
%int((A^2 - r^2)^2/r, r)
eq229 = D*S*C*Omega^2*(A^2 - B^2)/R/T/((A^4/4 - A^2*A^2 + A^4*log(A)) - (B^4/4 -...
A^2*B^2 + A^4*log(B)))*sigma;
simplify(Vpssolution(2)) - simplify(eq229)