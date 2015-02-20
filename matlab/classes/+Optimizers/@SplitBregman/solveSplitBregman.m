function u = solveSplitBregman(obj)
    % Grab local vars
    A = obj.A;
    At = obj.At;
    W = obj.W;
    Wt = obj.Wt;
    f = obj.f;
    u0 = obj.u0;
    Niter = obj.globalIter;
    cgIter = obj.cgIter;
    cgTol = obj.cgTol;
    lam = obj.lam;
    mu = obj.mu;
    
    % Initialize some local objects
    alpha = W(u0);
    v = alpha; %Lagrange multiplier. Copy construct OK since FrameletExpansion is not a handle class
    Atf = At(f);
    
    % Set up CGLS solver
    b = Atf + mu*Wt(alpha-v);
    cgSolver = Optimizers.CGLS(A,At,b,u0,mu,cgIter,cgTol);

    k = 0;
    while k<Niter     
        k=k+1;
        % Step 1: Approximate u(k+1) by solving (PtP+mu*I)u=Ptf0+muWt(alpha-v) 
        % using conjugate gradient.  
        u = cgSolver.solveCGLS();
        Wu = W(u);
        % Step 2: find alpha(k+1) by performing framelet soft thresholding
        alpha = v + Wu;
        alpha.softThreshold(lam/mu);
        % Step 3: update Lagrange multiplier
        %v = v + mu*(Wu-alpha);
        v = v + (Wu-alpha);
        % Step 4: update CGLS solver
        b = Atf + mu*Wt(alpha-v);
        clear cgSolver;
        cgSolver = Optimizers.CGLS(A,At,b,u,mu,cgIter,cgTol);
    end
    
  
end
    