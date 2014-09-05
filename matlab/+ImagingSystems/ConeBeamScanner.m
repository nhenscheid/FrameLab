classdef ConeBeamScanner
   properties (SetAccess = private)
       scanType = '3DCone'; %Support for other geometries later
       Na = 256;
       Nb = 256;
       Nv = 96;
       SO = single(54.1);%Source to isocenter
       OD = single(40.8);%Detector to isocenter
       dyDet; %Detector spacing
       dzDet; %Detector spacing
       sdPhi; 
       sdZ;
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
    end%Static methods
end