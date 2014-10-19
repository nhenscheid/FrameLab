classdef FrameletSystem
    properties(SetAccess=private)
        dim = 2;
        type = 'haar';
        level = 1;
        D;
        R; 
    end
    
    methods
        %***Constructor***%
        function obj = FrameletSystem(dim,type,level)
            if nargin>0
                %!!!Need to check validity of params
                obj.dim = dim;
                obj.type = type;
                obj.level = level;
                switch type
                    case 'haar'
                        [obj.D,obj.R] = obj.GenerateFrameletFilter(0);
                    case 'linear'
                        [obj.D,obj.R] = obj.GenerateFrameletFilter(1);
                    case 'cubic'
                        [obj.D,obj.R] = obj.GenerateFrameletFilter(3);
                    otherwise
                        %Throw some error
                end%Switch type
            end
        end%Constructor

        function alpha = forwardTransform(obj,u)
            %forwardTransform computes the framelet transform of the array
            %u.  The output is an object of type FrameletExpansion.  u is
            %a dim-dimensional Matlab array (should check that
            %dim(u)=obj.dim!!)
            switch obj.dim
                case 2
                    alpha = obj.FraDecMultiLevel2D(u,obj.D,obj.level);
                case 3
                    alpha = obj.FraDecMultiLevel3D(u,obj.D,obj.level);
                otherwise
                    %Throw some error
            end
            alpha = DataTypes.FrameletExpansion(obj.dim,obj.type,obj.level,alpha);
       
        end
        
        function u = adjointTransform(obj,alpha)
           %adjointTransform computes the adoint framelet transform of
           %alpha.  alpha is a cell array (FrameletExpansion.frameletArray)
           %The output is a Matlab array of dimension obj.dim 
           %(should check that dim(alpha) = obj.dim!!)
           switch obj.dim
                case 2
                    u = obj.FraRecMultiLevel2D(alpha,obj.R,obj.level);
                case 3
                    u = obj.FraRecMultiLevel3D(alpha,obj.R,obj.level);
                otherwise
                    %Throw some error
           end
        end
    end%Methods
    
    methods(Static,Access=private)
        [D,R]= GenerateFrameletFilter(frame);
        alpha = FraDecMultiLevel2D(A,D,L);
        alpha = FraDecMultiLevel3D(A,D,L);
        u = FraRecMultiLevel2D(alpha,R,L);
        u = FraRecMultiLevel3D(alpha,R,L);
    end
end