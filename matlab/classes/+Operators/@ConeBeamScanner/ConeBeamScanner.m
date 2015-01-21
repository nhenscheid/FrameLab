classdef ConeBeamScanner < handle
    %obj = ConeBeamScanner(type,na,nb,nv,zmax,rps,vtab,fps,phaseShift)
    % TO-DO (1/12/15)
    % Check for unbounded objects?
    properties (SetAccess = private)
        type = 'circle'; % circle or helix
        na = 256;
        nb = 256;
        nv = 128;
        SO = single(5.0);%Source to isocenter (cm)
        OD = single(5.0);%Detector to isocenter (cm)
        Ly = single(6.0); % Width of detector panel (cm)
        Lz = single(6.0); % Height of detector panel (cm)
        P = single(0.5); % Helix pitch (cm) 
        para = struct; %struct of parameters to pass to Gao's methods
        y_os = single(0.0);
        rps = 3; % Gantry rotations per second
        zmax = 2; % Scan length (cm)
        vtab = 2; % Table velocity (cm/s)
        fps = 30; % Frames per second
        nHelix = 1;
        phaseShift = 0; %
    end
    
    properties
        verbose = false; %To display diagnostic messages
        GPU = 1;  % Set to 0 to force using CPU
    end
    
    methods
        %***Constructor***%
        function obj = ConeBeamScanner(type,na,nb,nv,zmax,rps,vtab,fps,nHelix,phaseShift)
            %ConeBeamScanner(type,na,nb,nv,zmax,rps,vtab,fps)
            %!!!Should check inputs and throw errors if invalid
            if nargin < 4
                error('You must specify type, na, nb and nv');
            end
            if ~(strcmp(type,'circle')||strcmp(type,'helix')||strcmp(type,'multiHelix'))
                error('circle, helix and multiHelix are only valid types');
            end
            if nargin>0
                obj.type = type;
                obj.na = na;
                obj.nb = nb;
                obj.nv = nv;
            end
            if nargin>4 && strcmp(type,'circle')
                error('zmax, rps, vtab, fps are for helical scan only');
            end
            if nargin<5&& (strcmp(type,'helix')||strcmp(type,'multiHelix'))
                error('helix and multiHelix require more inputs')
            end
            if nargin>4 && (strcmp(type,'helix')||strcmp(type,'multiHelix'))
                if strcmp(type,'helix')
                    obj.nv = floor(fps*zmax/vtab)+1;
                elseif strcmp(type,'multiHelix')
                    obj.nv = nHelix*(floor(fps*zmax/vtab)+1);
                end
                obj.zmax = zmax;
                if nargin>5
                    obj.rps = rps;
                    if nargin>6
                        obj.vtab = vtab;
                        if nargin>7
                            obj.fps = fps;
                            if nargin>8
                                obj.nHelix = nHelix;
                                obj.phaseShift = phaseShift;
                            end
                        end
                    end
                end
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