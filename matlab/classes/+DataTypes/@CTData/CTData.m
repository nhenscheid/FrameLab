classdef CTData
    properties (SetAccess = private)
        scanner; %options are 'cone' and 'fan'
        dataArray; %array of size [na nb nv] 
        dataArrayNorm;
        L; %Object size 
    end
    
    
    methods
        %***Constructor***%
        function obj = CTData(scanner,dataArray,dataArrayNorm,L)
            if(nargin>0)
                obj.scanner = scanner;
                obj.dataArray = dataArray;
                obj.dataArrayNorm = dataArrayNorm;
                obj.L = L;
            end
        end % Constructor
        
        function alpha = frameletTransform(obj,sys)
            %frameletTransform computes the framelet transform of the
            %present object's data array using the FrameletSystem system.
            alpha = sys.forwardTransform(obj.dataArray);
        end%frameletTransform    
        
        function plotData3D(this,n,m,offset)
            if(strcmp(this.type,'helix')||strcmp(this.type,'multiHelix'))
                plotphantom3D(this.dataArray,n,m,offset);
            end
        end
        
        function plotData(this)
            if(length(size(this.dataArray))==2)
                imshow(this.dataArray,[]);
            end
        end
        
        function c = plus(a,b)
            %add an array to a CTData object
            if(size(a.dataArray)~=size(b))
                error('Array dimension mismatch!')
            end
            c = DataTypes.CTData(a.scanType,a.dataArray+b,a.para,a.L);
        end
        
        %John's equation
        D = applyJohn(obj)
    end
    
end