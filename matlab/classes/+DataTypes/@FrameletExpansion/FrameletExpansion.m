classdef FrameletExpansion
%classdef FrameletExpansion < handle & matlab.mixin.Copyable
   properties (SetAccess = private)
      dim = 3;
      frameletSystem;
      frameletArray;
   end
   
   methods
       %***Constructor***%
       function obj = FrameletExpansion(dim,frameletSystem,frameletArray)
           if nargin>0
               obj.dim = dim;
               obj.frameletSystem = frameletSystem;
               if nargin<3
                   %User did not supply frameletArray, so make an empty one
                   level = frameletSystem.level;
                   switch type
                       case 'haar'
                           obj.frameletArray = obj.makeCellArray(dim,2,level);
                       case 'linear'
                           obj.frameletArray = obj.makeCellArray(dim,3,level);
                       case 'cubic'
                           obj.frameletArray = obj.makeCellArray(dim,5,level);
                       otherwise
                           %Throw some error  
                   end%Switch type
               end
               if nargin==3
                   %User supplied a frameletArray, so use it
                   obj.frameletArray = frameletArray;
               end
           end
       end%Constructor 
       
       function u = adjointFrameletTransform(obj)
           %adjointFrameletTransform computes u = W^T\alpha where \alpha is
           %the framelet expansion data in this object and W^T is specified
           %by the FrameletTransform object sys.
           u = obj.frameletSystem.adjointTransform(obj.frameletArray);
       end
       
       plot(this)
       
       % Operators
       c = plus(a,b);
       c = minus(a,b);
       c = mtimes(a,b);    
       
       function plotData(this,n,m,offset)
            %currently specialized to 3D.
            plotphantom3D(this.dataArray,n,m,offset);
       end
       
   end%Methods 
   
   methods(Static)
       function X = makeCellArray(dim,size,level)
           X = cell(1,level);
           switch dim
               case 2
                   for i = 1:level
                   X{i} = cell(size,size);
                   end                   
               case 3
                   for i = 1:level
                   X{i} = cell(size,size,size);
                   end
               
               otherwise
                   error('Something is wrong.');
           end%Switch dim
       end%makeCellArray
       
   end%Static methods
   
    
    
    
end