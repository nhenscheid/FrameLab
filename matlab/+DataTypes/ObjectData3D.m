classdef ObjectData3D
    properties (SetAccess = private)
        meshType = 'cartesian'; %later: perhaps triangulated etc?
        Nx = 256; 
        Ny = 256; 
        Nz = 256; 
        Lx = 10;
        Ly = 10;
        Lz = 10;
        dx;
        dy;
        dz;
        dataArray;
    end %Properties
    
    methods
        %***Constructor***%
        function obj = ObjectData3D(array,Lx,Ly,Lz)
            %!!!Should check inputs and throw errors if invalid.
            if nargin >0
                arrSize = size(array);
                obj.dataArray = array;
                obj.Nx = arrSize(1);
                obj.Ny = arrSize(2);
                obj.Nz = arrSize(3);
                obj.Lx = Lx;
                obj.Ly = Ly;
                obj.Lz = Lz;
                obj.dx = Lx/arrSize(1);
                obj.dy = Ly/arrSize(2);
                obj.dz = Lz/arrSize(3); 
            else
                obj.dx = obj.Lx/obj.Nx;
                obj.dy = obj.Ly/obj.Ny;
                obj.dz = obj.Lz/obj.Nz;
            end
        end%Constructor
        
        function coord = getCoords(obj,i,j,k)
            %getCoords(i,j,k) computes the (x,y,z) coordinates of the
            %(i,j,k)th voxel.
            coord = [0,0,0];
            coord(1) = -obj.Lx/2+obj.dx*(i-0.5);
            coord(2) = -obj.Ly/2+obj.dy*(j-0.5);
            coord(3) = -obj.Lz/2+obj.dz*(k-0.5);
        end%getCoords
        
        function alpha = frameletTransform(obj,sys)
            %frameletTransform computes the framelet transform of the
            %present object's data array using the FrameletSystem system.
            alpha = sys.forwardTransform(obj.dataArray);
        end%frameletTransform    
    end %Methods
end%Classdef
        