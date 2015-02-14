function [x,y,f1] = ivpTest(N,lam)
    h = 1/(N-1);
    f = @(x)cos(x);
    
    D = (diag(-1*ones(N-1,1),-1)+diag(ones(N-1,1),1))/(2*h);
    D(1,1) = -1/h;
    D(1,2) = 1/h;
    D(end,end-1) = -1/h;
    D(end,end) = 1/h;
    DtD = D'*D;
    
    RtR = zeros(N);
    RtR(1,1) = 1;
    %RtR(2,2) = 1;
    %RtR(end-1,end-1) = 1;
    RtR(end,end) = 1;
    
    rank(RtR+lam*DtD)

    x = (0:h:1)';
    f1 = f(x);
    %rhs = [1;f(h);zeros(N-4,1);f(1-h);f(1)]+lam*D'*f1;
    rhs = [1;zeros(N-2,1);f(1)]+lam*D'*f1;
    y = (RtR+lam*DtD)\rhs;
end



%function RtR(y)
    


%end