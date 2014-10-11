classdef filter 
    properties (SetAccess = private)
        a;
    end
    
    methods 
        function obj = filter(a)
            obj.a = a;
        end
        
        function objDat = apply(obj,u)
            objDat = copy(u);
            objDat.updateDataArray(obj.a*u.dataArray);
        end
        
        
        
    end
    
    
end