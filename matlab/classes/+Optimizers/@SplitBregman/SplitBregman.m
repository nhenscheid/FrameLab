classdef SplitBregman
   % The SplitBregman class solves problems of the type 
   % min_x (1/2)||Au-f||_2^2 + lam ||Wu||_1 where W is a tight frame and A
   % is a linear operator.
   
   properties (SetAccess = private)
      globalIter=1;
      cgIter=10;
      cgTol=1e-10;
      A;
      At;
      W;
      Wt;
      u0;
      f;
      lam;
      mu;
   end
   
   methods
      %***Constructor***%
      function obj = SplitBregman(A,At,W,Wt,lam,mu,f,u0,globalIter,cgIter,cgTol)
          if nargin<8
              error('SplitBregman requires A,At,W,Wt,lam,f,and u0');
          end
          if nargin>=8
             obj.A = A;
             obj.At = At;
             obj.W = W;
             obj.Wt = Wt;
             obj.lam = lam;
             obj.mu = mu;
             obj.f = f;
             obj.u0 = u0;
             if nargin >=9
                 obj.globalIter = globalIter;
                 if nargin>=10
                     obj.cgIter = cgIter;
                     if nargin>=11
                         obj.cgTol = cgTol;
                     end
                 end
             end           
          end
      end %Constructor
       
      % Solve system using Split Bregman
      u = solveSplitBregman(obj)
   end
    
    
    
end