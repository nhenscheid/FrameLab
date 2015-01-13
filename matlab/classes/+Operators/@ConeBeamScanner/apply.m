function y=apply(this,object)
    this.setPara(object); %set the Gao parameter struct
    X0 = object.dataArray(:);
    if (gpuDeviceCount==0 || this.GPU == 0)
        this.para.GPU = uint32(0);
        this.checkInputs(X0);
        if(this.verbose)
            disp('Computing forward cone beam transform with CPU');
            this.para
        end
        [y,yNorm] = Ax_cone_mf_cpu(X0,this.para);
    elseif (gpuDeviceCount>0 && this.GPU ==1)
        this.para.GPU = uint32(1);
        this.checkInputs(X0);
        if(this.verbose)
            disp('Computing forward cone beam transform with GPU');
            this.para
        end
        y = Ax_cone_mf(X0,this.para);
    end
    if strcmp(this.type,'circle')
        y = DataTypes.CTData('circle',reshape(y,[this.na,this.nb,this.nv]),this.para,object.L,reshape(yNorm,[this.na,this.nb,this.nv]));
    elseif strcmp(this.type,'helix')
        y = DataTypes.CTData('helix',reshape(y,[this.na this.nb this.nv]),this.para,object.L,reshape(yNorm,[this.na,this.nb,this.nv]));
    elseif strcmp(this.type,'doubleHelix')
        size(y)
        this.na
        this.nb
        this.nv
        y1 = y(1:this.na*this.nb*this.nv/2);
        y1Norm = yNorm(1:this.na*this.nb*this.nv/2); 
        y2 = y(this.na*this.nb*this.nv/2+1:end);
        y2Norm = yNorm(this.na*this.nb*this.nv/2+1:end);
        size(y1)
        size(y2)
        A = zeros(this.na,this.nb,this.nv/2,2);
        ANorm = zeros(this.na,this.nb,this.nv/2,2);
        A(:,:,:,1) = reshape(y1,[this.na this.nb this.nv/2]);
        A(:,:,:,2) = reshape(y2,[this.na this.nb this.nv/2]);
        ANorm(:,:,:,1) = reshape(y1Norm,[this.na this.nb this.nv/2]);
        ANorm(:,:,:,2) = reshape(y2Norm,[this.na this.nb this.nv/2]);        
        y = DataTypes.CTData('doubleHelix',A,this.para,object.L,ANorm);
    end
end