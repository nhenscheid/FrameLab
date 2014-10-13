classdef CGLS
    %The CGLS class solves problems of the type Ax = b where A is a
    %symmetric positive definite linear operator.  The method is matrix
    %free - x and b can be matlab objects and A a linear operator mapping
    %objects of type x to type b.
    properties (SetAccess = private)
        globalIter;
        convergenceTol;
        lhs; %lhs is the operator - a function handle 
        rhs;
    end


    methods 
        %***Constructor***%
        function obj = CGLS(lhs,rhs,globalIter,convergenceTol)
            if nargin>0
                obj.lhs = lhs;
                obj.rhs = rhs;
                obj.globalIter = globalIter;
                obj.convergenceTol = convergenceTol;
            end
        end %constructor
        
        %Solve system using cgls
        y = solveCGLS(obj)
    end




end