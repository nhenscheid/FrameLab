% Testing the ObjectData class

disp('Creating empty ObjectData')
u = DataTypes.ObjectData()
disp('Creating empty 2D ObjectData')
u = DataTypes.ObjectData(2)
disp('Creating empty 3D ObjectData')
u = DataTypes.ObjectData(3)
disp('Creating 2D ObjectData with zeros')
u = DataTypes.ObjectData(2,zeros(256))
disp('Creating 3D ObjectData with zeros')
u = DataTypes.ObjectData(3,zeros(256,256,256))