function y = changeScanner(ctData,objectData,dir)
    % Not sure how else to do this, so we're just going to make a function.
    % y will be a CTData object with new scanner settings
    % For now, we're going to change back and forth between a 2 helix scan
    % and a 1 helix scan.
    
    scanner1 = ctData.scanner;
    na = scanner1.na;
    nb = scanner1.nb;
    nv = scanner1.nv/scanner1.nHelix;
    zmax = scanner1.zmax;
    rps = scanner1.rps;
    vtab = scanner1.vtab;
    fps = scanner1.fps;
    dphi = 2*pi*rps/fps;
    
    if strcmp(dir,'down')
        nHelix = 2;
        phaseShift = dphi*(0:nHelix-1);
        scanner2 = Operators.ConeBeamScanner('multiHelix',na,nb,nv,zmax,rps,vtab,fps,nHelix,phaseShift);
        scanner2.setPara(objectData);
        dataArray = single(zeros(na,nb,nv,2));
        dataArrayNorm = single(zeros(na,nb,nv,2));
        dataArray(:,:,:,1) = ctData.dataArray(:,:,:,1);
        dataArray(:,:,:,2) = ctData.dataArray(:,:,:,3);
        dataArrayNorm(:,:,:,1) = ctData.dataArrayNorm(:,:,:,1);
        dataArrayNorm(:,:,:,2) = ctData.dataArrayNorm(:,:,:,3);
        y = DataTypes.CTData(scanner2,dataArray,dataArrayNorm,ctData.L);
    elseif strcmp(dir,'up')
        nHelix = 3;
        phaseShift = dphi*(0:nHelix-1); 
        scanner2 = Operators.ConeBeamScanner('multiHelix',na,nb,nv,zmax,rps,vtab,fps,nHelix,phaseShift);
        scanner2.setPara(objectData);
        dataArray = single(zeros(na,nb,nv,3));
        dataArrayNorm = single(zeros(na,nb,nv,3));
        dataArray(:,:,:,1) = ctData.dataArray(:,:,:,1);
        dataArray(:,:,:,3) = ctData.dataArray(:,:,:,2);
        dataArrayNorm(:,:,:,1) = ctData.dataArrayNorm(:,:,:,1);
        dataArrayNorm(:,:,:,3) = ctData.dataArrayNorm(:,:,:,2);
        y = DataTypes.CTData(scanner2,dataArray,dataArrayNorm,ctData.L);
    else
        error('dir can only be up or down');
    end
end