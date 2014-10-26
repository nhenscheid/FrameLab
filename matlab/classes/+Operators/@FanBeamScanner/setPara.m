function setPara(this,object)
    %setPara builds the parameter struct to pass to Hao Gao's CT methods
    
    %!!!NOTE: Rescaling is to force dx = dy = 1, as per Gao 2012
    scale = single(object.Dx(1));
    
    %%%%Scanner parameters***%
    this.para.SO = this.SO/scale; %Source-isocenter spacing
    this.para.OD = this.OD/scale; %Isocenter-detector plane spacing
    this.para.dy_det = this.dy_det/scale; %Detector pixel spacing
    this.para.y_os = single(0.0);
   
    %***Object parameters***%
    N = size(object.dataArray);
    this.para.scale = scale; %Object xy scaling factor to force dx=dy=1
    this.para.nx = uint32(N(1)); %Object array x-dim
    this.para.ny = uint32(N(2)); %Object array y-dim
    this.para.nt = uint32(1); %Number of time steps (set to 0 for now)
    
    %***Scan variables (derived from scanner params)
    this.para.sd_phi = single(2*pi/this.nv*(0:this.nv-1));
    this.para.y_det=single(((-this.nd/2:this.nd/2-1)+0.5)*this.dy_det+this.y_os)/scale;
    this.para.cos_phi = cos(this.para.sd_phi);
    this.para.sin_phi = sin(this.para.sd_phi);
    this.para.cos_det = [];
    this.para.sin_det = [];
    angle_det=atan2(this.para.y_det,this.para.SO+this.para.OD);
    this.para.cos_det=cos(angle_det);
    this.para.sin_det=sin(angle_det);
    
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
    this.para.nv = uint32(this.nv);
    this.para.tmp_size = uint32(tmp_size);
    %this.para.nv_block = uint32(4); %dont need this for 2d?
    this.para.version = uint32(1); %
    

end%setPara