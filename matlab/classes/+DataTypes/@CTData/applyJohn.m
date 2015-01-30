function D = applyJohn(obj,smoothed)
    %D=applyJohn(CTData,smoothed) 
    %Apply John's equation for a multiHelix scan
    %Smoothed = 1 if we apply to a smoothed version, 0 else
    % Extract some scanner parameters 
    cbct = obj.scanner;
    nHelix = cbct.nHelix;
    if(nHelix<2)
        error('Johns equation only applies on multiHelix data!');
    end
    para = cbct.para;
    scale = para.scale;
    nv = para.Nv/nHelix;
    na = cbct.na;
    nb = cbct.nb;
    SO = cbct.SO;
    OD = cbct.OD;
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
    R = SO+OD;
    dzeta = -phaseShift(2)*h %This assumes that each phase shift is the same?
    dphi = 2*pi*cbct.rps/cbct.fps
    da = scale*cbct.para.dy_det
    db = scale*cbct.para.dz_det
    
    % Initialize derivative arrays
%     gtb = single(zeros(size(g)));
%     gaz = gtb;
%     gbz = gtb;
%     gb = gtb;
%     gbb = gtb;
%     gab = gtb;

    % gtb
    gtb = (g(:,3:end,3:end,:)+g(:,1:end-2,1:end-2,:)-g(:,1:end-2,3:end,:)-g(:,3:end,1:end-2,:))/(4*db*dphi);
    % gaz, gbz
    if(nHelix<3)
        gaz = (g(3:end,:,:,2)+g(1:end-2,:,:,1)-g(1:end-2,:,:,2)-g(3:end,:,:,1))/(2*da*dzeta);
        gbz = (g(:,3:end,:,2)+g(:,1:end-2,:,1)-g(:,1:end-2,:,2)-g(:,3:end,:,1))/(2*db*dzeta);
    else
        gaz = (g(3:end,:,:,3:end)+g(1:end-2,:,:,1:end-2)-g(1:end-2,:,:,3:end)-g(3:end,:,:,1:end-2))/(4*da*dzeta);
        gbz = (g(:,3:end,:,3:end)+g(:,1:end-2,:,1:end-2)-g(:,1:end-2,:,3:end)-g(:,3:end,:,1:end-2))/(4*db*dzeta);
    end
    % gb
    gb = (g(:,3:end,:,:)-g(:,1:end-2,:,:))/(2*db);
    % gbb
    gbb = (g(:,3:end,:,:)+g(:,1:end-2,:,:)-2*g(:,2:end-1,:,:))/(db*db);
    % gab
    gab = (g(3:end,3:end,:,:)+g(1:end-2,1:end-2,:,:)-g(3:end,1:end-2,:,:)-g(1:end-2,3:end,:,:))/(4*da*db);
    % a
    [a,b,NULL1,NULL2] = ndgrid(y_det,z_det,zeros(nv,1),zeros(nHelix,1));
    % b 

    size(gtb)
    size(gaz)
    size(gbz)
    size(gb)
    size(gbb)
    size(gab)
    size(a)
    size(b)
    if(nHelix<3)
        D = R*gtb(2:end-1,:,:,1) - R*SO*gaz(:,2:end-1,2:end-1) -...
            R*h*gbz(2:end-1,:,2:end-1) +...
            2.*a(2:end-1,2:end-1,2:end-1,1).*gb(2:end-1,:,2:end-1,1) +...
            a(2:end-1,2:end-1,2:end-1,1).*b(2:end-1,2:end-1,2:end-1,1).*gbb(2:end-1,:,2:end-1,1) +...
            (a(2:end-1,2:end-1,2:end-1,1).^2+R^2).*gab(:,:,2:end-1,1);
    else
        D = R*gtb(2:end-1,:,:,2:end-1) - R*SO*gaz(:,2:end-1,2:end-1,:) -...
            R*h*gbz(2:end-1,:,2:end-1,:) +...
            2.*a(2:end-1,2:end-1,2:end-1,2:end-1).*gb(2:end-1,:,2:end-1,2:end-1) +...
            a(2:end-1,2:end-1,2:end-1,2:end-1).*b(2:end-1,2:end-1,2:end-1,2:end-1).*gbb(2:end-1,:,2:end-1,2:end-1) +...
            (a(2:end-1,2:end-1,2:end-1,2:end-1).^2+R^2).*gab(:,:,2:end-1,2:end-1);
    end
    
    D = DataTypes.CTData(cbct,obj.dataArray,single(D),obj.L); 

end