classdef CGLS
    
    properties (SetAccess = private)
        globalIter;
        convergenceTol;
        lhs; 
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
        
        function 
        
        
        
    end




end