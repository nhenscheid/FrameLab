classdef ConeBeamData
    properties (SetAccess = private)
        scanType; %spiral, circular, ?
        Na;
        Nb;
        Nv;
        SO = single(54.1);
        OD = single(40.8);
        dyDet;
        dzDet;
        sdPhi;
        sdZ;
        dataArray;
    end
    
    
    methods
        %***Constructor***%
        function obj = ConeBeamData()
            
            
            
        end
        
        
    end
    
    
    
    
    
    
end