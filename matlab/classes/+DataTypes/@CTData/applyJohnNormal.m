function D = applyJohnNormal(obj,smoothed)
    %D=applyJohnNormal(CTData,smoothed) 
    %Apply DtD for John's equation, for a multiHelix scan
    %Smoothed = 1 if we apply to a smoothed version, 0 else
    % Extract some scanner parameters 
    cbct = obj.scanner;
    nz = double(cbct.nHelix);
    if(nz<2)
        error('Johns equation only applies on multiHelix data!');
    end
    para = cbct.para;
    scale = para.scale;
    nv = double(para.Nv/nz);
    na = double(cbct.na);
    nb = double(cbct.nb);
    SO = double(cbct.SO);
    OD = double(cbct.OD);
    cos_phi = para.cos_phi;
    sin_phi = para.sin_phi;
    sd_z = scale*para.sd_z;
    h = cbct.vtab/(2*pi*cbct.rps); 
    y_det = scale*(cbct.para.y_det);
    z_det = scale*(cbct.para.z_det);
    phaseShift = cbct.phaseShift;
    
    %Warning: specialized to 2 helix data!
    %g1 = obj.dataArrayNorm(:,:,:,1);
    %g2 = obj.dataArrayNorm(:,:,:,2);
    g = obj.dataArrayNorm;

    if(smoothed)
        framelet = Transforms.FrameletSystem(3,'linear',1);
        disp('Computing framelet expansion of f')
        alpha = framelet.forwardTransform(g);
        g = alpha.frameletArray{1}{1,1};
        %alpha1 = framelet.forwardTransform(g1);
        %alpha2 = framelet.forwardTransform(g2);
        %g1 = alpha1.frameletArray{1}{1,1};
        %g2 = alpha2.frameletArray{1}{1,1};
    end

    % Compute approximate John's equation 
    disp('approximating Johns Equation')
    R = double(SO+OD);
    dzeta = double(-phaseShift(2)*h); %This assumes that each phase shift is the same?
    dtheta = double(2*pi*cbct.rps/cbct.fps);
    da = double(scale*cbct.para.dy_det);
    db = double(scale*cbct.para.dz_det);
    
    % Upgrade g to a double matrix (Matlab can only do sparse doubles) 
    g = double(g);
    
    % Build a and b 
    [a,b,NULL1,NULL2] = ndgrid(y_det,z_det,zeros(nv,1),zeros(nz,1));
    
    % Construct the 1D finite difference matrices
    Da1 = spdiags([-ones(na,1),zeros(na,1),ones(na,1)],0:2,na-2,na)/(2*da);
    Db1 = spdiags([-ones(nb,1),zeros(nb,1),ones(nb,1)],0:2,nb-2,nb)/(2*db);
    Dbb1 = spdiags([ones(nb,1),-2*ones(nb,1),ones(nb,1)],0:2,nb-2,nb)/(db^2);
    Dt1 = spdiags([-ones(nv,1),zeros(nv,1),ones(nv,1)],0:2,nv-2,nv)/(2*dtheta);
    if nz>2
        Dz1 = spdiags([-ones(nz,1),zeros(nz,1),ones(nz,1)],0:2,nz-2,nz)/(2*dzeta);
    elseif nz==2
        Dz1 = spdiags([-ones(nz,1),ones(nz,1)],0:1,nz-1,nz)/dzeta;
    end
    % Construct the 4D finite difference matrices
    Dtb = kron(kron(speye(nz),Dt1),kron(Db1,speye(na)));
    Daz = kron(kron(Dz1,speye(nv)),kron(speye(nb),Da1));
    Dbz = kron(kron(Dz1,speye(nv)),kron(Db1,speye(na)));
    Db = kron(kron(speye(nz),speye(nv)),kron(Db1,speye(na)));
    Dbb = kron(kron(speye(nz),speye(nv)),kron(Dbb1,speye(na)));
    Dab = kron(kron(speye(nz),speye(nv)),kron(Db1,Da1));
    
    % Build some trim matrices
    aTrim = speye(na);aTrim = aTrim(2:end-1,:);
    bTrim = speye(nb);bTrim = bTrim(2:end-1,:);
    tTrim = speye(nv);tTrim = tTrim(2:end-1,:);
    zTrim = speye(nz);
        if nz>2
            zTrim = zTrim(2:end-1,:);
            azTrim = kron(kron(zTrim,speye(nv-2)),kron(speye(nb-2),aTrim));
            btTrim = kron(kron(speye(nz-2),tTrim),kron(bTrim,speye(na-2)));
            atTrim = kron(kron(speye(nz-2),tTrim),kron(speye(nb-2),aTrim));
            atzTrim = kron(kron(zTrim,tTrim),kron(speye(nb-2),aTrim));
            tzTrim = kron(kron(zTrim,tTrim),kron(speye(nb-2),speye(na-2)));
        else
            zTrim = zTrim(1,:);
            azTrim = kron(kron(zTrim,speye(nv-2)),kron(speye(nb-2),aTrim));
            btTrim = kron(kron(speye(nz-1),tTrim),kron(bTrim,speye(na-2)));
            atTrim = kron(kron(speye(nz-1),tTrim),kron(speye(nb-2),aTrim));
            atzTrim = kron(kron(zTrim,tTrim),kron(speye(nb-2),aTrim));
            tzTrim = kron(kron(zTrim,tTrim),kron(speye(nb-2),speye(na-2)));
        end
        
  
    nbInt = na*(nb-2)*nv*nz;
    nabInt = (na-2)*(nb-2)*nv*nz; 
    a1 = double(spdiags(reshape(a(:,2:end-1,:,:),[],1),0,nbInt,nbInt));
    b1 = double(spdiags(reshape(b(:,2:end-1,:,:),[],1),0,nbInt,nbInt));
    a2 = double(spdiags(reshape(a(2:end-1,2:end-1,:,:).^2 + R^2,[],1),0,nabInt,nabInt));
    
    
    Dfull = R*azTrim*Dtb - R*SO*btTrim*Daz - R*h*atTrim*Dbz + 2*atzTrim*a1*Db + atzTrim*a1*b1*Dbb + tzTrim*a2*Dab;
%     disp('Size of Dfull is ')
%     size(Dfull)
%     
%     whos
%     % gtb
%     ByteSize(Dtb)
%     ByteSize(Daz)
%     ByteSize(Dbz)
%     ByteSize(Db)
%     ByteSize(Dbb)
%     ByteSize(Dab)
%     ByteSize(Dfull)
%     
    
    gtb = reshape(Dtb*g(:),[na,nb-2,nv-2,nz]);
    
    % gaz, gbz
    if(nz<3)     
        gaz = reshape(Daz*g(:),[na-2,nb,nv,nz-1]);
        gbz = reshape(Dbz*g(:),[na,nb-2,nv,nz-1]);
    else
        gaz = reshape(Daz*g(:),[na-2,nb,nv,nz-2]);
        gbz = reshape(Dbz*g(:),[na,nb-2,nv,nz-2]);
    end
    % gb
    gb = reshape(Db*g(:),[na,nb-2,nv,nz]);
    % gbb
    gbb = reshape(Dbb*g(:),[na,nb-2,nv,nz]);
    % gab
    gab = reshape(Dab*g(:),[na-2,nb-2,nv,nz]);
    % a and b
    
    
    if(nz<3)
        %D = R*gtb(2:end-1,:,:,1) - R*SO*gaz(:,2:end-1,2:end-1,1) -...
        %    R*h*gbz(2:end-1,:,2:end-1,1) +...
        %    2.*a(2:end-1,2:end-1,2:end-1,1).*gb(2:end-1,:,2:end-1,1) +...
        %    a(2:end-1,2:end-1,2:end-1,1).*b(2:end-1,2:end-1,2:end-1,1).*gbb(2:end-1,:,2:end-1,1) +...
        %    (a(2:end-1,2:end-1,2:end-1,1).^2+R^2).*gab(:,:,2:end-1,1);
        Dnew = reshape(Dfull*g(:),[na-2,nb-2,nv-2,1]);
        %disp('norm(D-Dnew)')
        %norm(D(:)-Dnew(:))/norm(D(:))
    else
        D = R*gtb(2:end-1,:,:,2:end-1) - R*SO*gaz(:,2:end-1,2:end-1,:) -...
            R*h*gbz(2:end-1,:,2:end-1,:) +...
            2.*a(2:end-1,2:end-1,2:end-1,2:end-1).*gb(2:end-1,:,2:end-1,2:end-1) +...
            a(2:end-1,2:end-1,2:end-1,2:end-1).*b(2:end-1,2:end-1,2:end-1,2:end-1).*gbb(2:end-1,:,2:end-1,2:end-1) +...
            (a(2:end-1,2:end-1,2:end-1,2:end-1).^2+R^2).*gab(:,:,2:end-1,2:end-1);
    end
    
   D = reshape((Dfull')*Dnew(:),[na,nb,nv,nz]);

end