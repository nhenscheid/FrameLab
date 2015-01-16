function y = coneBeamSphere(t,a,b,zeta,h,Rs,Rd,normalized)
% y = coneBeamSphere(t,a,b,z,h,Rs,Rd,normalized)
% Helical cone beam transform of unit ball
% t ~ gantry angle
% (a,b) ~ detector coords
% zeta ~ axial shift (scalar)
% h = P/2pi, P ~ helix pitch
% Rs,Rd ~ source to iso, detector to iso distances
% normalized ~ bool, should we normalize by ||xs-xd||?
    R = Rs+Rd;
    if(~normalized)
        y = 2*sqrt(((Rs*R-b.*(h*t+zeta)).^2-(a.^2+b.^2+R^2).*(Rs^2+(h*t+zeta+1).*(h*t+zeta-1)))./(a.^2+b.^2+R^2));
    else
        y = 2*sqrt((Rs*R-b.*(h*t+zeta)).^2-(a.^2+b.^2+R^2).*(Rs^2+(h*t+zeta+1).*(h*t+zeta-1)))./(a.^2+b.^2+R^2);
    end
        y(imag(y)~=0)=0;
end