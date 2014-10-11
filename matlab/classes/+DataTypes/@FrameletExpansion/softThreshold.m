function softThreshold(obj,tau)
    %softThreshold computes T_zeta(alpha) where alpha is this
    %object's framelet array.
    switch obj.dim
        case 2
            temp = shrinkFramelet2D(obj.frameletArray,tau);
        case 3
            temp = shrinkFramelet3D(obj.frameletArray,tau); 
        otherwise
            %Throw an error?
    end
        obj.frameletArray = temp;
end

function R=shrinkFramelet2D(A,tau)
    %===================================================
    %
    %  Execute the shrinkage of Framelet coefficients
    %  stored in array A with threshhold tau. The low
    %  "frequencies" coefficients are keep without any
    %  changes.
    %  
    %  Parameters:
    %       A: input array of Framelet coeeficients
    %       tau: threshhold
    %       R: output array of threshholded coefficients
    %
    %   Author: Jerome Gilles
    %   Institution: UCLA - Math Department
    %   email: jegilles@math.ucla.edu
    %====================================================
    L=length(A);
    for l=1:L
       [NH,NW]=size(A{l});
       for nh=1:NH
          for nw=1:NW
              if ((nw == 1) && (nh == 1))
                 R{l}{nh,nw}=A{l}{nh,nw};
              else
                 R{l}{nh,nw}=sign(A{l}{nh,nw}).*max(zeros(size(A{l}{nh,nw})),abs(A{l}{nh,nw})-tau);
              end
          end
       end
    end
end%shrinkFramelet2D

function R=shrinkFramelet3D(A,tau)
    %===================================================
    %
    %  Execute the shrinkage of Framelet coefficients
    %  stored in array A with threshhold tau. The low
    %  "frequencies" coefficients are keep without any
    %  changes.
    %  
    %  Parameters:
    %       A: input array of Framelet coeeficients
    %       tau: threshhold
    %       R: output array of threshholded coefficients
    %
    %   Author: Jerome Gilles
    %   Institution: UCLA - Math Department
    %   email: jegilles@math.ucla.edu
    %====================================================
    disp(['shrinking 3D framelet array with tau = ',num2str(tau)])
    L=length(A);
    for l=1:L
       [NH,NW]=size(A{l});
       for ii=1:NH
          for jj=1:NH
              for kk=1:NH
                  if ((jj == 1) && (ii == 1) && (kk==1))
                      R{l}{ii,jj,kk}=A{l}{ii,jj,kk};
                  else
                      R{l}{ii,jj,kk}=sign(A{l}{ii,jj,kk}).*max(zeros(size(A{l}{ii,jj,kk})),abs(A{l}{ii,jj,kk})-tau);
                  end
              end%for kk=1:NH
          end%forjj=1:NH
       end%for ii=1:NH
    end%for i=1:L
end%shrinkFramelet3D