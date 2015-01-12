classdef CTData
    properties (SetAccess = private)
        scanType; %options are 'cone' and 'fan'
        para; %Parameter struct
        dataArray; %array of size [na nb nv] 
        L; %Object size 
    end
    
    
    methods
        %***Constructor***%
        function obj = CTData(scanType,dataArray,para,L)
            if nargin>1
               dataArray = single(dataArray); 
            end
            if(nargin<4)
               scanType = 'fan'; 
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
        
        function alpha = frameletTransform(obj,sys)
            %frameletTransform computes the framelet transform of the
            %present object's data array using the FrameletSystem system.
            alpha = sys.forwardTransform(obj.dataArray);
        end%frameletTransform    
        
        function plotData3D(this,n,m,offset)
            validatestring(this.scanType,{'cone'});
            plotphantom3D(this.dataArray,n,m,offset);
        end
        
        function plotData(this)
            validatestring(this.scanType,{'fan'});
            imshow(this.dataArray,[]);
        end
        
        function c = plus(a,b)
            %add an array to a CTData object
            if(size(a.dataArray)~=size(b))
                error('Array dimension mismatch!')
            end
            c = DataTypes.CTData(a.scanType,a.dataArray+b,a.para,a.L);
        end
    end
    
end