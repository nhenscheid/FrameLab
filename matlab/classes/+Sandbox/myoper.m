classdef myoper
   properties
      B = zeros(10); 
      alpha = 1.0
   end
   
   
   methods
       %***Constructor***%
       function obj = myoper(B,alpha)
           if nargin>0
               obj.B = B;
               obj.alpha = alpha;
           end
       end%Constructor
       
       function C = doOper(this,A)
           %This function takes the matrix A and outputs the matrix 
           %C = alpha*B*A
           C = this.alpha*this.B*A;
       end
       
   end
    
    
end