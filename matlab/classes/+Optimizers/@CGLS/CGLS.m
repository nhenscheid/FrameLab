classdef CGLS
    % The CGLS class solves problems of the type Ax = b using regularized
    % least squares, i.e. x*:= (AtA+lam*I)\At(b)
    % The method is matrix free - x and b can be matlab objects and 
    % A a linear operator mapping objects of type x to type b.
    % CGLS uses regularized least squares - i.e. we solve 
    properties (SetAccess = private)
        globalIter;
        convergenceTol;
        A; % The forward operator (a function handle)
        At; % The adjoint operator (a function handle)
        b; % Rhs
        u0; % Initial solution guess 
        lam; % Regularization parameter
        uExact;
    end


    methods 
        %***Constructor***%
        function obj = CGLS(A,At,b,u0,lam,globalIter,convergenceTol,uExact)
            if nargin<5
               error('CGLS requires A, At, b, u0 and lam'); 
            end
            if nargin>=5
                obj.A = A;
                obj.b = b;
                obj.u0 = u0;
                obj.At = At;
                obj.lam = lam;
                if nargin==6
                    obj.globalIter = globalIter;
                    obj.convergenceTol = 1e-14;
                elseif nargin==7
                    obj.globalIter = globalIter;
                    obj.convergenceTol = convergenceTol;
                elseif nargin==8
                    obj.uExact = uExact;
                else
                    obj.globalIter = 500; %Prob too big
                    obj.convergenceTol = 1e-14; %Prob too small 
                end
            end
           
        end %constructor
        
        %Solve system using cgls
        y = solveCGLS(obj)
    end




end