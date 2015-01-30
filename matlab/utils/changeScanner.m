function y = changeScanner(ctData,objectData,dir)
    % Not sure how else to do this, so we're just going to make a function.
    % y will be a CTData object with new scanner settings
    % For now, we're going to change back and forth between a 2 helix scan
    % and a 1 helix scan.
    
    scanner1 = ctData.scanner;
    na = scanner1.na;
    nb = scanner1.nb;
    nv = scanner1.nv;
    zmax = scanner1.zmax;
    rps = scanner1.rps;
    vtab = scanner1.vtab;
    fps = scanner1.fps;
    dphi = 2*pi*rps/fps;
    
    if strcmp(dir,'down')
        nHelix = 1;
        phaseShift = 0;
        scanner2 = Operators.ConeBeamScanner('multiHelix',na,nb,nv,zmax,rps,vtab,fps,nHelix,phaseShift)
        scanner2.setPara(objectData);
        y = DataTypes.CTData(scanner2,ctData.dataArray(:,:,:,1),ctData.dataArrayNorm(:,:,:,1),ctData.L);
    elseif strcmp(dir,'up')
        nHelix = 2;
        phaseShift = dphi*(0:nHelix-1); 
        scanner2 = Operators.ConeBeamScanner('multiHelix',na,nb,nv,zmax,rps,vtab,fps,nHelix,phaseShift)
        scanner2.setPara(objectData);
        dataArray = single(zeros(na,nb,nv,2));
        dataArrayNorm = single(zeros(na,nb,nv,2));
        dataArray(:,:,:,1) = ctData.dataArray;
        dataArrayNorm(:,:,:,1) = ctData.dataArrayNorm;
        y = DataTypes.CTData(scanner2,dataArray,dataArrayNorm,ctData.L);
    else
        error('dir can only be up or down');
    end
end