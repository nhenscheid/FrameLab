classdef ConeBeamData
    properties (SetAccess = private)
        scanType; %spiral, circular, ?
        para; %Parameter struct
        dataArray; %array of size [na nb nv] 
        L; %Object size 
    end
    
    
    methods
        %***Constructor***%
        function obj = ConeBeamData(dataArray,para,L,scanType)
            if(nargin<4)
               scanType = 'default'; 
            end
            obj.checkInputs(dataArray,para,scanType);
            if(nargin>0)
                obj.para = para;
                obj.dataArray = dataArray;
                obj.L = L;
                if(nargin==4)
                    obj.scanType = scanType;
                end            
            end
        end % Constructor
        
        function plotData(this,n,m,offset)
           plotphantom3D(this.dataArray,n,m,offset);
        end
    end
    
end