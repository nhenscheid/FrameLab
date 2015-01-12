function y = solveCGLS(obj)
% TO-DO: 
% - Pre-conditioning? 
% - Input checks/warnings/etc
% - error update/verbose mode
    lam = obj.lam;
    A_loc = @(x)(obj.At(obj.A(x))+lam*x);
    %rhs = obj.At(obj.b);
    rk = A_loc(obj.u0)-obj.b;
    pk = -rk;
    uk = obj.u0;
    gammaold = rk.objdot(rk); %object rk must have dot product method
    icg = 0;
    cgiter = obj.globalIter;
    cgtol = obj.convergenceTol;
    % Main conjugate gradient loop
    while icg<cgiter
        icg = icg+1;
        qk = A_loc(pk);
        alphak = gammaold/(pk.objdot(qk));
        uk = uk + alphak*pk;
        rk = rk + alphak*qk;
        gammanew = rk.objdot(rk)
        betak = gammanew/gammaold;
        pk = -rk + betak*pk;
        gammaold = gammanew;
        if(gammaold<cgtol)
            icg = cgiter;
        end
    end
    y = uk;
end