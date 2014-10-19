classdef FrameletExpansion < handle
   properties (SetAccess = private)
      dim = 3;
      frameletType = 'haar';
      nLevel = 1;
      frameletArray;
   end
   
   methods
       %***Constructor***%
       function obj = FrameletExpansion(dim,type,level,frameletArray)
           if nargin>0
               obj.dim = dim;
               obj.frameletType = type;
               obj.nLevel = level;
               if nargin==3
                   %User did not supply frameletArray, so make an empty one
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
               if nargin==4
                   %User supplied a frameletArray, so use it
                   obj.frameletArray = frameletArray;
               end
           end
       end%Constructor 
       
       function u = adjointFrameletTransform(obj,sys)
           %adjointFrameletTransform computes u = W^T\alpha where \alpha is
           %the framelet expansion data in this object and W^T is specified
           %by the FrameletTransform object sys.
           u = sys.adjointTransform(obj.frameletArray);
       end
       
       c = plus(a,b);
           
       
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
                   %Throw some error
           end%Switch dim
       end%makeCellArray
       
   end%Static methods
   
    
    
    
end