function plotphantom3D(p,n,m,offset)
% Plot a grid of 32 slices of the 3D phantom p.

%figure
ha=tight_subplot(n,m);
for iplot=1:n*m
    axes(ha(iplot));imshow(p(:,:,iplot+offset),[])
end