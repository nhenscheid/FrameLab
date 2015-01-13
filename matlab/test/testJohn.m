% Testing John's equation, helix parametrization

clear all;
N = 256;
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]);
%u0 = DataTypes.ObjectData(3,single(gaussian3D(N)),[2,2,2]);
nd = 128;
nv = 256;
phase = pi/32;

cbct = Operators.ConeBeamScanner('doubleHelix',nd,nd,[],1,4,1,nv,phase);
cbct.verbose = true;
cbct.GPU = 1;

disp('Computing double helix scan of u0')
f = cbct.apply(u0);

% Extract some scanner parameters 
nv = cbct.nv/2;
na = cbct.na;
nb = cbct.nb;
SO = cbct.SO;
scale = cbct.para.scale;
cos_phi = cbct.para.cos_phi;
sin_phi = cbct.para.sin_phi;
sd_z = scale*cbct.para.sd_z;

% Compute detector positions
OD = cbct.OD;
y_det = scale*(cbct.para.y_det);
z_det = scale*(cbct.para.z_det);

% Verify that f is correct
h = cbct.vtab/(2*pi*cbct.rps); 
% l = @(t,a1,a2)2*sqrt(-(a1^2*(24+h^2*t^2)+4*(-25+(2*a2+5*h*t)*(3*a2+5*h*t)))/(100+a1^2+a2^2));
% yexact = zeros(nd,nd,nv);
% 
% disp('finding exact solution for helix 1')
% for iv = 1:nv
%     for i = 1:nd
%         for j = 1:nd
%             t = cbct.para.sd_phi(iv);
%             a1 = y_det(i);
%             a2 = z_det(j);
%             if imag(l(t,a1,a2))==0
%                 yexact(i,j,iv) = l(t,a1,a2);
%             end
%         end
%     end
% end
% 
% y = f.dataArray(:,:,:,1);
% disp(sprintf('%s%d\n','L2 error for first helix is ',sqrt(sum(abs(yexact(:)-y(:)).^2))/sqrt(sum(abs(yexact(:)).^2))))
% disp(sprintf('%s%d\n','sup error for first helix is ',max(abs(yexact(:)-y(:)))/max(abs(yexact(:)))))
% 
% l2 = @(t,a1,a2)sqrt(-384-(h*(pi+4*t))^2+((-200+a2*h*(pi+4*t))^2)/(100+a1^2+a2^2))/2;
% yexact2 = zeros(nd,nd,nv);
% 
% disp('finding exact solution for helix 2')
% for iv = 1:nv
%     for i = 1:nd
%         for j = 1:nd
%             t = cbct.para.sd_phi(iv);
%             a1 = y_det(i);
%             a2 = z_det(j);
%             if imag(l(t,a1,a2))==0
%                 yexact2(i,j,iv) = l2(t,a1,a2);
%             end
%         end
%     end
% end
% 
% y = f.dataArray(:,:,:,2);
% disp(sprintf('%s%d\n','L2 error for second helix is ',sqrt(sum(abs(yexact2(:)-y(:)).^2))/sqrt(sum(abs(yexact2(:)).^2))))
% disp(sprintf('%s%d\n','sup error for second helix is ',max(abs(yexact2(:)-y(:)))/max(abs(yexact2(:)))))
% 
% % ***Compute normalized transform***
% % Source and detector positions
% xSource =@(k)[SO*cos_phi(k);SO*sin_phi(k);sd_z(k)]; 
% xDet=@(i,j,k)[-OD*cos_phi(k)-y_det(i)*sin_phi(k);-OD*sin_phi(k)+y_det(i)*cos_phi(k);z_det(j)+sd_z(k)];

%scan1 = f.dataArrayNorm(:,:,:,1);
%scan2 = f.dataArrayNorm(:,:,:,2);

%g1 = scan1;
%g2 = scan2;
% Form scaled transform arrays

% disp('scaling arrays');
% for i=1:na
%     for j=1:nb
%         for k=1:nv
%             g1(i,j,k) = scan1(i,j,k)/norm(xSource(k)-xDet(i,j,k));
%             g2(i,j,k) = scan2(i,j,k)/norm(xSource(k+nv)-xDet(i,j,k+nv));
%         end
%     end
% end

g1 = f.dataArrayNorm(:,:,:,1);
g2 = f.dataArrayNorm(:,:,:,2);

% lNorm = @(t,a,b,z)(2*sqrt( (-24-(h*t+z)^2)/(100+a^2+b^2) + (-50+b*(h*t+z)^2)^2/(100+a^2+b^2)^2));
% 
% yexactNorm = zeros(nd,nd,nv);
% 
% disp('finding exact solution for normalized ray transform')
% for iv = 1:nv
%     for i = 1:nd
%         for j = 1:nd
%             t = cbct.para.sd_phi(iv);
%             a1 = y_det(i);
%             a2 = z_det(j);
%             if imag(l(t,a1,a2))==0
%                 yexactNorm(i,j,iv) = lNorm(t,a1,a2,0);
%             end
%         end
%     end
% end
% 
% disp(sprintf('%s%d\n','L2 error for second helix is ',sqrt(sum(abs(yexactNorm(:)-g1(:)).^2))/sqrt(sum(abs(yexactNorm(:)).^2))))
% disp(sprintf('%s%d\n','sup error for second helix is ',max(abs(yexactNorm(:)-g1(:)))/max(abs(yexactNorm(:)))))



framelet = Transforms.FrameletSystem(3,'linear',1);

disp('Computing framelet expansion of u')
alpha1 = framelet.forwardTransform(g1);
alpha2 = framelet.forwardTransform(g2);

g1 = alpha1.frameletArray{1}{1,1};
g2 = alpha2.frameletArray{1}{1,1};


% Compute approximate John's equation 
disp('approximating Johns Equation')
R = SO+OD;
%dzeta = 0.03125;  % NOTE: THIS IS SPECIFIC TO CURRENT SCAN SETTINGS, NOTABLY THE PI/4 PHASE SHIFT
dzeta = (cbct.para.sd_z(2)-cbct.para.sd_z(1))*cbct.para.scale;
dphi = 2*pi*cbct.rps/cbct.fps;
da = scale*cbct.para.dy_det;
db = scale*cbct.para.dz_det;


% gtb
gtb = (g1(:,3:end,3:end)+g1(:,1:end-2,1:end-2)-g1(:,1:end-2,3:end)-g1(:,3:end,1:end-2))/(4*db*dphi);
% gaz
gaz = (g2(3:end,:,:)+g1(1:end-2,:,:)-g2(1:end-2,:,:)-g1(3:end,:,:))/(4*da*dzeta);
% gbz
gbz = (g2(:,3:end,:)+g1(:,1:end-2,:)-g2(:,1:end-2,:)-g1(:,3:end,:))/(4*db*dzeta);
% gb
gb = (g1(:,3:end,:)-g1(:,1:end-2,:))/(2*db);
% gaa
gaa = (g1(3:end,:,:)+g1(1:end-2,:,:)-2*g1(2:end-1,:,:))/(da*da);
% gab
gab = (g1(3:end,3:end,:)+g1(1:end-2,1:end-2,:)-g1(3:end,1:end-2,:)-g1(1:end-2,3:end,:))/(4*da*db);
% a
[a,b,NULL] = ndgrid(y_det,z_det,zeros(nv,1));
% b 


D = R*gtb(2:end-1,:,:) - R*SO*gaz(:,2:end-1,2:end-1) -...
    R*h*gbz(2:end-1,:,2:end-1) +...
    2.*a(2:end-1,2:end-1,2:end-1).*gb(2:end-1,:,2:end-1) +...
    a(2:end-1,2:end-1,2:end-1).*b(2:end-1,2:end-1,2:end-1).*gaa(:,2:end-1,2:end-1) +...
    (a(2:end-1,2:end-1,2:end-1).^2+R^2).*gab(:,:,2:end-1);







