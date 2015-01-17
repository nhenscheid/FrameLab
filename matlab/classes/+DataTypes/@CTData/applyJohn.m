% Implementation of John's equation for helical scans
function D = applyJohn(obj)
    % Extract some scanner parameters 
    cbct = obj.scanner;
    nHelix = cbct.nHelix;
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
    g1 = obj.dataArrayNorm(:,:,:,1);
    g2 = obj.dataArrayNorm(:,:,:,2);

    framelet = Transforms.FrameletSystem(3,'linear',2);
    disp('Computing framelet expansion of f')
    alpha1 = framelet.forwardTransform(g1);
    alpha2 = framelet.forwardTransform(g2);

    g1 = alpha1.frameletArray{2}{1,1};
    g2 = alpha2.frameletArray{2}{1,1};

    % Compute approximate John's equation 
    disp('approximating Johns Equation')
    R = SO+OD;
    dzeta = -phaseShift(2)*h
    dphi = 2*pi*cbct.rps/cbct.fps
    da = scale*cbct.para.dy_det
    db = scale*cbct.para.dz_det

    % gtb
    gtb = (g1(:,3:end,3:end)+g1(:,1:end-2,1:end-2)-g1(:,1:end-2,3:end)-g1(:,3:end,1:end-2))/(4*db*dphi);
    % gaz
    gaz = (g2(3:end,:,:)+g1(1:end-2,:,:)-g2(1:end-2,:,:)-g1(3:end,:,:))/(4*da*dzeta);
    % gbz
    gbz = (g2(:,3:end,:)+g1(:,1:end-2,:)-g2(:,1:end-2,:)-g1(:,3:end,:))/(4*db*dzeta);
    % gb
    gb = (g1(:,3:end,:)-g1(:,1:end-2,:))/(2*db);
    % gbb
    gbb = (g1(:,3:end,:)+g1(:,1:end-2,:)-2*g1(:,2:end-1,:))/(db*db);
    % gab
    gab = (g1(3:end,3:end,:)+g1(1:end-2,1:end-2,:)-g1(3:end,1:end-2,:)-g1(1:end-2,3:end,:))/(4*da*db);
    % a
    [a,b,NULL] = ndgrid(y_det,z_det,zeros(nv,1));
    % b 

    D = R*gtb(2:end-1,:,:) - R*SO*gaz(:,2:end-1,2:end-1) -...
        R*h*gbz(2:end-1,:,2:end-1) +...
        2.*a(2:end-1,2:end-1,2:end-1).*gb(2:end-1,:,2:end-1) +...
        a(2:end-1,2:end-1,2:end-1).*b(2:end-1,2:end-1,2:end-1).*gbb(2:end-1,:,2:end-1) +...
        (a(2:end-1,2:end-1,2:end-1).^2+R^2).*gab(:,:,2:end-1);

end