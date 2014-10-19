classdef ObjectData < handle & matlab.mixin.Copyable %not sure why handle
    properties (SetAccess = private)
        dim;
        meshtype = 'cartesian'; %only option for now
        N;  %array dimensions [Nx,Ny] or [Nx,Ny,Nz]
        L;  %object size [Lx,Ly] or [Lx,Ly,Lz]
        Dx; %voxel size 
        dataArray;
    end
    
    methods 
        %***Constructor***%
        function obj = ObjectData(dim,array,L)
            %!!!INPUT CHECK!!!%
            if nargin>=1
                obj.dim = dim;
            end
            if nargin==1
                switch dim
                    case 2
                        obj.N = [128,128];
                        obj.L = [10,10];
                        obj.Dx = [10,10]./[128,128];
                        obj.dataArray = zeros(128);
                    case 3
                        obj.N = [128,128,128];
                        obj.L = [10,10,10];
                        obj.Dx = [10,10,10]./[128,128,128];
                        obj.dataArray = zeros(128,128,128);
                end%switch dim
            elseif nargin==2
                obj.dataArray = array;
                obj.N = size(array);
                switch dim
                    case 2
                        obj.L = [10,10];
                        obj.Dx = [10,10]./size(array);
                    case 3
                        obj.L = [10,10,10];
                        obj.Dx = [10,10,10]./size(array);
                end
            elseif nargin==3
                obj.dim = dim;
                obj.dataArray = array;
                obj.N = size(array);
                obj.L = L;
                obj.Dx = L./size(array); %only makes sense for cart grid
            end
        end%Constructor
        
        function getCoords(obj,X)
            %Some function to return the coordinates of a voxel
        end
        
        function alpha = frameletTransform(obj,sys)
            %frameletTransform computes the framelet transform of the
            %present object's data array using the FrameletSystem system.
            alpha = sys.forwardTransform(obj.dataArray);
        end%frameletTransform    
        
        function adjointFrameletTransform(obj,sys,alpha)
            %adjointFrameletTransform replaces this object's data array
            %with W^T\alpha, where W is specified by the FrameletTransform
            %object sys and alpha is an object of FrameletExpansion type.
            obj.dataArray = alpha.adjointFrameletTransform(sys);
        end%adjointFrameletTransform
        
        function f = applyOperator(obj,oper)
            %applyOperator takes this object and applies the operator oper
            %to it.  This is very abstract - oper eats objects of type
            %ObjectData and spits out *some* kind of object. For example,
            %it might be an object of type ObjectData, or of type
            %ConeBeamData, etc.
            f = oper.apply(obj);
        end
        
        function updateDataArray(obj,A)
            %updateDataArray is a property setter method for obj.dataArray
            obj.dataArray = A; 
        end
        
        function plotData(this,n,m,offset)
            %currently specialized to 3D.
            plotphantom3D(this.dataArray,n,m,offset);
        end
        
        function obj = uminus(obj)
            %unary minus, returns the same object with negative data
            %entries
            obj.dataArray = -obj.dataArray;
        end
        
        function obj = mtimes(alpha,obj)
            %multiply this object's dataArray by a matrix
            obj.dataArray = alpha*obj.dataArray;
        end
        
        function c = plus(a,b)
            %add two objectDatas and store the result in the first's dataArray
            if((a.dim~=b.dim)||(max(a.L~=b.L)))
                error('Either dimensions or sizes are wrong!')
            end
            c = DataTypes.ObjectData(a.dim,a.dataArray+b.dataArray,a.L);
        end
        
        function y = objdot(obj,obj2)
            %dot product of two objects of type ObjectData
            y = dot(obj.dataArray(:),obj2.dataArray(:));
        end
        
        
        
        
    end%methods

end%classdef