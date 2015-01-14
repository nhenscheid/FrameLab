%Test Adjoint

N = 128;
na = 128;
nb = 128;
nv = 64;

disp('testing adjoint for circular scan')

cbct = Operators.ConeBeamScanner('circle',na,nb,nv);
cbct.GPU = 1;

% Compute random 3D array and apply cone beam transform 
x = DataTypes.ObjectData(3,rand(N,N,N),[2,2,2]);
Ax = cbct.apply(x);

% Compute random CT data set and apply adjoint 
y = DataTypes.CTData('cone',rand(na,nb,nv),cbct.para,[2,2,2]);
Aty = cbct.applyAdjoint(y);

% Test <Ax,y> = <x,Aty>
% (Should write a dot method for these data types)
abs(dot(Ax.dataArray(:),y.dataArray(:))-dot(x.dataArray(:),Aty.dataArray(:)))/abs(dot(Ax.dataArray(:),y.dataArray(:)))


disp('testing adjoint for helical scan')

cbct = Operators.ConeBeamScanner('helix',na,nb,[],0.5,3,0.25,64);
cbct.GPU = 1;

% Compute random 3D array and apply cone beam transform 
x = DataTypes.ObjectData(3,rand(N,N,N),[2,2,2]);
Ax = cbct.apply(x);

% Compute random CT data set and apply adjoint 
y = DataTypes.CTData('cone',rand(na,nb,cbct.nv),cbct.para,[2,2,2]);
Aty = cbct.applyAdjoint(y);

% Test <Ax,y> = <x,Aty>
% (Should write a dot method for these data types)
abs(dot(Ax.dataArray(:),y.dataArray(:))-dot(x.dataArray(:),Aty.dataArray(:)))/abs(dot(Ax.dataArray(:),y.dataArray(:)))
abs(dot(Ax.dataArray(:),y.dataArray(:))-dot(x.dataArray(:),Aty.dataArray(:)))


disp('testing adjoint for helical scan with phantom')

cbct = Operators.ConeBeamScanner('helix',na,nb,[],0.5,3,0.25,64);
cbct.GPU = 1;

% Compute random 3D array and apply cone beam transform 
x = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);
Ax = cbct.apply(x);

% Compute random CT data set and apply adjoint 
y = DataTypes.CTData('cone',rand(na,nb,cbct.nv),cbct.para,[2,2,2]);
Aty = cbct.applyAdjoint(y);

% Test <Ax,y> = <x,Aty>
% (Should write a dot method for these data types)
abs(dot(Ax.dataArray(:),y.dataArray(:))-dot(x.dataArray(:),Aty.dataArray(:)))/abs(dot(Ax.dataArray(:),y.dataArray(:)))
abs(dot(Ax.dataArray(:),y.dataArray(:))-dot(x.dataArray(:),Aty.dataArray(:)))

