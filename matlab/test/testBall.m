% eTest the cone beam scanner with a unit ball
clear all;clc
N = 128;
u0 = DataTypes.ObjectData(3,single(ball3d(N)),[2,2,2]);
nd = 16;

cbct = Operators.ConeBeamScanner(nd,nd,1);

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

disp(sprintf('%s%d\n','Max error is ',max(abs(yexact(:)-y.dataArray(:)))))
disp(sprintf('%s%d\n','L2 error is ',sqrt(sum(abs(yexact(:)-y.dataArray(:)).^2))))