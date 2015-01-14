% Test framelet expansion
%u = DataTypes.ObjectData(3,single(phantom3d(256)),[2,2,2]);
u = DataTypes.ObjectData(3,rand(10,10,10),[2,2,2])
disp('Creating 3D 1-level Haar wavelet transform system')
framelet = Transforms.FrameletSystem(3,'haar',1);

disp('Computing framelet expansion of u')
alpha = u.frameletTransform(framelet);

disp('Computing adjoint framelet expansion of alpha = Wu')

WTWu = alpha.adjointFrameletTransform(framelet);