A = rand(3,2)
rank(A)
b = rand(3,1)
mlsol = A\b


cgiter = 10;
cgtol = 1e-12;

u = zeros(2,1);
rk = b;
pk = rk;
P = @(x)A*x;
Pt = @(x)A'*x;
mu = 0; %for regularization
gammaold = norm(rk(:))^2;
icg = 0;

while icg < cgiter
            icg=icg+1;
            qk=P(Pt(pk))+mu*pk;
            alphak = gammaold/dot(pk(:),qk(:));
            u = u + alphak*pk;
            rk = rk + alphak*qk;
            gammanew = norm(rk(:))^2;
            betak = gammanew/gammaold;
            pk = -rk + betak*pk;
            gammaold = gammanew;
            if (gammaold<cgtol)
                icg=cgiter;
            end
end