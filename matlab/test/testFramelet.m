% Test framelet expansion
type = input('What type? (haar, linear or cubic) ','s')

if(~(strcmp(type,'haar')||strcmp(type,'linear')||strcmp(type,'cubic')))
    error('Incorrect type - only haar, linear or cubic allowed')
end
disp('Testing 2D Framelet Expansion')
u2D = DataTypes.ObjectData(2,phantom(128),[2,2]);
disp('Creating 2D 1-level wavelet transform system')
framelet = Transforms.FrameletSystem(2,type,1);
disp('Computing framelet expansion');
alpha2D = framelet.forwardTransform(u2D);
alpha2D.plot();

Wtalpha2D = alpha2D.adjointFrameletTransform();
disp('Checking for perfect reconstruction: ')
norm(Wtalpha2D(:) - u2D.dataArray(:))/norm(u2D.dataArray(:))

u3D = DataTypes.ObjectData(3,phantom3d(128),[2,2,2]);
%u = DataTypes.ObjectData(3,rand(10,10,10),[2,2,2])
disp('Creating 3D 1-level Haar wavelet transform system')
framelet = Transforms.FrameletSystem(3,type,1);

disp('Computing framelet expansion of u')
alpha3D = framelet.forwardTransform(u3D);

disp('Computing adjoint framelet expansion of alpha = Wu')

Wtalpha3D = alpha3D.adjointFrameletTransform();
disp('Checking for perfect reconstruction: ')
norm(Wtalpha3D(:)- u3D.dataArray(:))/norm(u3D.dataArray(:))