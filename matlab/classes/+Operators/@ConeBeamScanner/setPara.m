function setPara(this,object)
    %setPara builds the parameter struct to pass to Hao Gao's CT methods
    
    %!!!NOTE: Rescaling is to force dx = dy = 1, as per Gao 2012
    scale = single(object.Dx(1));
    
    %%%%Scanner parameters***%
    this.para.SO = this.SO/scale; %Source-isocenter spacing
    this.para.OD = this.OD/scale; %Isocenter-detector plane spacing
    this.para.dy_det = this.dy_det/scale; %Detector pixel spacing
    this.para.dz_det = this.dz_det/scale; %Detector pixel spacing
    this.para.y_os = this.y_os/this.para.dy_det;
   
    %***Object parameters***%
    N = size(object.dataArray);
    this.para.scale = scale; %Object xy scaling factor to force dx=dy=1
    this.para.nx = uint32(N(1)); %Object array x-dim
    this.para.ny = uint32(N(2)); %Object array y-dim
    this.para.nz = uint32(N(3)); %Object array z-dim
    this.para.dz = single(object.Dx(3))/scale; %Object dz
    this.para.nt = uint32(1); %Number of time steps (set to 0 for now)
    
    %***Scan variables (derived from scanner params)
    this.para.sd_phi = single(2*pi/this.nv*(0:this.nv-1));
    this.para.sd_z = single(zeros(1,this.nv))/scale;%for helical scan
    this.para.y_det=single(((-this.na/2:this.na/2-1)+0.5)*this.dy_det+this.y_os)/scale;
    this.para.z_det=single(((-this.nb/2:this.nb/2-1)+0.5)*this.dz_det)/scale;
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

end%setPara