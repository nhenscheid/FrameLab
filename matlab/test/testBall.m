% Test the fan and cone beam scanners with a unit disc/ball
%% Fan beam geometry
N = 10000;
u0 = DataTypes.ObjectData(2,single(unitBall(N,2)),[2,2]);
nd = 10;
nv = 1;
fbct = Operators.FanBeamScanner(nd,nv)
y = fbct.apply(u0);
y = y.dataArray;
y_det = fbct.para.y_det*fbct.para.scale;

l = @(a1)4*sqrt(-6+625/(100+a1^2));

yexact = zeros(nd,1);
for i=1:nd
    if imag(l(y_det(i)))==0
        yexact(i) = l(y_det(i));
    end
end
disp(sprintf('%s%d\n','Max error is ',max(abs(y-yexact(:)))/max(abs(yexact))))
disp(sprintf('%s%d\n','L2 error is ',norm(y-yexact)/norm(yexact)))
%% Circular source path
clear all;clc
N = 256;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
nd = 128;
nv = 1;

disp('Testing circular scan');
cbct = Operators.ConeBeamScanner('circle',nd,nd,nv);

cbct.verbose = true;
cbct.GPU = 1;

y = cbct.apply(u0);

y_det = cbct.para.y_det*cbct.para.scale;
z_det = y_det;

% Exact solution (z = 0 slice) 
l=@(a1,a2)4*sqrt(-6+625/(100+a1^2+a2^2));

yexact = zeros(nd,nd);

for i = 1:nd
    for j = 1:nd
        a1 = y_det(i);
        a2 = z_det(j);
        if imag(l(a1,a2))==0
            yexact(i,j) = l(y_det(i),z_det(j));
        end
    end
end

for i = 1:nv
    disp('For view ');
    i
    y1 = y.dataArray(:,:,i);
    disp(sprintf('%s%d\n','Max error is ',max(abs(yexact(:)-y1(:)))))
    disp(sprintf('%s%d\n','L2 error is ',sqrt(sum(abs(yexact(:)-y1(:)).^2))/sqrt(sum(abs(yexact(:))))))
end

%% Helical 
disp('Testing helical scan');
cbct = Operators.ConeBeamScanner('helix',nd,nd,[],2,4,2,64);
cbct.verbose = true;
cbct.GPU = 1;

y = cbct.apply(u0);
y_det = cbct.para.y_det*cbct.para.scale;
z_zet = y_det;

h = cbct.vtab/(2*pi*cbct.rps); 
l = @(t,a1,a2)2*sqrt(-(a1^2*(24+h^2*t^2)+4*(-25+(2*a2+5*h*t)*(3*a2+5*h*t)))/(100+a1^2+a2^2));
nv = cbct.para.Nv;
yexact = zeros(nd,nd,nv);

for iv = 1:nv
    for i = 1:nd
        for j = 1:nd
            t = cbct.para.sd_phi(iv);
            a1 = y_det(i);
            a2 = z_det(j);
            if imag(l(t,a1,a2))==0
                yexact(i,j,iv) = l(t,a1,a2);
            end
        end
    end
end


disp(sprintf('%s%d\n','L2 error is ',sqrt(sum(abs(yexact(:)-y.dataArray(:)).^2))/sqrt(sum(abs(yexact(:))))))

