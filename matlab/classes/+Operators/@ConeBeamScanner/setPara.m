function setPara(this,object)
    %setPara builds the parameter struct to pass to Hao Gao's CT methods
    
    %!!!NOTE: Rescaling is to force dx = dy = 1, as per Gao 2012
    scale = single(object.Dx(1));
    
    %%%%Scanner parameters***%
    this.para.SO = this.SO/scale; %Source-isocenter spacing
    this.para.OD = this.OD/scale; %Isocenter-detector plane spacing
    dy_det = this.Ly/this.na; %Detector pixel spacing
    dz_det = this.Lz/this.nb; %Detector pixel spacing
    this.para.dy_det = dy_det/scale;
    this.para.dz_det = dz_det/scale;
    this.para.y_os = this.y_os/this.para.dy_det;
   
    %***Object parameters***%
    N = size(object.dataArray);
    this.para.scale = scale; %Object xy scaling factor to force dx=dy=1
    this.para.nx = uint32(N(1)); %Object array x-dim
    this.para.ny = uint32(N(2)); %Object array y-dim
    this.para.nz = uint32(N(3)); %Object array z-dim
    this.para.dz = single(object.Dx(3))/scale; %Object dz
    this.para.nt = uint32(1); %Number of time steps (set to 1 for now)
    
    %***Scan variables (derived from scanner params)
    if strcmp(this.type,'circle')
        this.para.sd_phi = single(2*pi/this.nv*(0:this.nv-1));
        this.para.sd_z = single(zeros(1,this.nv))/scale;
    elseif strcmp(this.type,'helix')
        phiMax = this.rps*this.zmax*2*pi/this.vtab;
        dphi = 2*pi*this.rps/this.fps;
        this.para.sd_phi = single(-phiMax/2:dphi:phiMax/2);
        h = this.vtab/(2*pi*this.rps);
        %this.para.sd_z = single(h*this.para.sd_phi-this.zmax/2)/scale;
        this.para.sd_z = single(h*this.para.sd_phi)/scale;
    elseif strcmp(this.type,'multiHelix')
        phiMax = this.rps*this.zmax*2*pi/this.vtab;
        dphi = 2*pi*this.rps/this.fps;
        nHelix = this.nHelix;
        this.para.sd_phi = [];
        for i=1:nHelix
           this.para.sd_phi = [this.para.sd_phi,single((-phiMax/2:dphi:phiMax/2)+this.phaseShift(i))];
        end
        %this.para.sd_phi = [single(0:dphi:phiMax),single((0:dphi:phiMax)+this.phaseShift)]; % second helix is phase shifted
        h = this.vtab/(2*pi*this.rps);
        %this.para.sd_z = [single(h*(0:dphi:phiMax)-this.zmax/2)/scale,single(h*(0:dphi:phiMax)-this.zmax/2)/scale]; % same z-coords for both
        this.para.sd_z = [];
        for i=1:nHelix
            this.para.sd_z = [this.para.sd_z,single(h*(-phiMax/2:dphi:phiMax/2))/scale];
        end
    end
    
    this.para.y_det=single(((-this.na/2:this.na/2-1)+0.5)*dy_det+this.y_os)/scale;
    this.para.z_det=single(((-this.nb/2:this.nb/2-1)+0.5)*dz_det)/scale;
    this.para.cos_phi = cos(this.para.sd_phi);
    this.para.sin_phi = sin(this.para.sd_phi);
    
    
    %***Misc things***%
    nt = 1;
    Id_v=cell(1,this.para.nt);
    Id_v{1}=uint32(0:this.nv-1);
    id_X=[];
    Nv=zeros(1,this.para.nt);
    for i=1:nt
        id_X=[id_X (i-1)*ones(1,numel(Id_v{i}))];
        Nv(i)=numel(Id_v{i});
    end
    tmp_size=max(Nv);
    id_Y=zeros(tmp_size,this.para.nt);
    for i=1:nt
        id_Y(1:Nv(i),i)=Id_v{i};
    end
    
    this.para.id_X = uint32(id_X);
    this.para.id_Y = uint32(id_Y);
    this.para.Nv = uint32(Nv);
    
    
    this.para.tmp_size = uint32(tmp_size);
    this.para.nv_block = uint32(4);
    this.para.version = uint32(1);  % 1 = Gao, 0 = Siddon

end%setPara