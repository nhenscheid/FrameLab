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
    
    
    methods
        %***Constructor***%
        function obj = ConeBeamScanner(Na,Nb,Nv)
            %!!!Should check inputs and throw errors if invalid
            if nargin>0
                obj.Na = Na;
                obj.Nb = Nb;
                obj.Nv = Nv;
            end
        end%Constructor
    end
    
    methods (Static = true)
        %***Geometry Plotter***%
        function plotGeometry()
            %Do nothing 
        end
        
        function y=doScan(this,object)
            this.setPara(object); %set the Gao parameter struct
            X0 = object.dataArray(:);
            this.checkInputs(X0);
            disp('computing forward cone beam transform!');
            size(X0)
            this.para
            y = Ax_cone_mf(X0,this.para);
        end
        

    end%Static methods
end