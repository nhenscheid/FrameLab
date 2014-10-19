classdef ConeBeamScanner < handle
    properties (SetAccess = private)
        scanType = '3DCone'; %Support for other geometries later
        na = 256;
        nb = 256;
        nv = 128;
        SO = single(54.1);%Source to isocenter
        OD = single(40.8);%Detector to isocenter
        dy_det =single(0.15/2); %Detector spacing  (!!!NEED DEFAULT VALUE!!!)
        dz_det =single(0.15/2); %Detector spacing  (!!!NEED DEFAULT VALUE!!!)
        para = struct; %struct of parameters to pass to Gao's methods
        y_os = single(0.0);
    end
    
    properties
        verbose = false; %To display diagnostic messages
    end
    
    
    methods
        %***Constructor***%
        function obj = ConeBeamScanner(na,nb,nv)
            %!!!Should check inputs and throw errors if invalid
            if nargin>0
                obj.na = na;
                obj.nb = nb;
                obj.nv = nv;
            end
        end%Constructor
        
        %***Forward scan***%
        y = apply(this,object)
        
        %***Adjoint operator***%
        y = applyAdjoint(this,object)
    end %Methods 
    
    methods (Static = true)
        %***Geometry Plotter***%
        function plotGeometry()
            %Do nothing 
        end
    end%Static methods
end