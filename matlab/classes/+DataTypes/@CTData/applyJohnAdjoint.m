% Implementation of John's equation for helical scans
function D = applyJohnAdjoint(obj,smoothed)
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

    if(smoothed)
        framelet = Transforms.FrameletSystem(3,'linear',1);
        disp('Computing framelet expansion of f')
        alpha1 = framelet.forwardTransform(g1);
        alpha2 = framelet.forwardTransform(g2);
        g1 = alpha1.frameletArray{1}{1,1};
        g2 = alpha2.frameletArray{1}{1,1};
    end
   

    % Compute approximate John's equation 
    disp('approximating Johns Equation')
    R = SO+OD;
    dzeta = -phaseShift(2)*h
    dphi = 2*pi*cbct.rps/cbct.fps
    da = scale*cbct.para.dy_det
    db = scale*cbct.para.dz_det

    % a and b
    [a,b,NULL] = ndgrid(y_det,z_det,zeros(nv,1));
    ag1 = a.*g1;
    abg1 = a.*b.*g1;
    a2g1 = (a.^2).*g1;
    % gtb
    g_tb = (g1(:,3:end,3:end)+g1(:,1:end-2,1:end-2)-g1(:,1:end-2,3:end)-g1(:,3:end,1:end-2))/(4*db*dphi);
    % gaz
    g_az = (g2(3:end,:,:)+g1(1:end-2,:,:)-g2(1:end-2,:,:)-g1(3:end,:,:))/(4*da*dzeta);
    % gbz
    g_bz = (g2(:,3:end,:)+g1(:,1:end-2,:)-g2(:,1:end-2,:)-g1(:,3:end,:))/(4*db*dzeta);
    % gb
    ag_b = (ag1(:,3:end,:)-ag1(:,1:end-2,:))/(2*db);
    % gbb
    abg_bb = (abg1(:,3:end,:)+abg1(:,1:end-2,:)-2*abg1(:,2:end-1,:))/(db*db);
    % gab
    g_ab = (g1(3:end,3:end,:)+g1(1:end-2,1:end-2,:)-g1(3:end,1:end-2,:)-g1(1:end-2,3:end,:))/(4*da*db);
    % (a^2g)ab
    a2g_ab = (a2g1(3:end,3:end,:)+a2g1(1:end-2,1:end-2,:)-a2g1(3:end,1:end-2,:)-a2g1(1:end-2,3:end,:))/(4*da*db);
    
    

    D = R*g_tb(2:end-1,:,:) - R*SO*g_az(:,2:end-1,2:end-1) -...
        R*h*g_bz(2:end-1,:,2:end-1) -...
        2.*ag_b(2:end-1,:,2:end-1) +...
        abg_bb(2:end-1,:,2:end-1) +...
        R^2*g_ab(:,:,2:end-1) +...
        a2g_ab(:,:,2:end-1);

end