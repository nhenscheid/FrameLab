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
        [y,yNorm] = Ax_cone_mf(X0,this.para);
    end
    if strcmp(this.type,'circle')
        y = DataTypes.CTData(this,reshape(y,[this.na,this.nb,this.nv]),reshape(yNorm,[this.na,this.nb,this.nv]),object.L);
    elseif strcmp(this.type,'helix')
        y = DataTypes.CTData(this,reshape(y,[this.na this.nb this.nv]),reshape(yNorm,[this.na,this.nb,this.nv]),object.L);
    elseif strcmp(this.type,'multiHelix')
        na = this.na;
        nb = this.nb;
        nv = this.nv/this.nHelix;
        nHelix = this.nHelix;
        
        A = single(zeros(na,nb,nv,nHelix));
        ANorm = single(zeros(na,nb,nv,nHelix));
        
        for i=1:nHelix
            j1 = (i-1)*na*nb*nv+1;
            j2 = na*nb*i*nv;
            A(:,:,:,i) = reshape(y(j1:j2),[na,nb,nv]);  
            ANorm(:,:,:,i) = reshape(yNorm(j1:j2),[na,nb,nv]); 
        end
%         y1 = y(1:this.na*this.nb*this.nv/2);
%         y1Norm = yNorm(1:this.na*this.nb*this.nv/2); 
%         y2 = y(this.na*this.nb*this.nv/2+1:end);
%         y2Norm = yNorm(this.na*this.nb*this.nv/2+1:end);
%         size(y1)
%         size(y2)
%         A = zeros(this.na,this.nb,this.nv/2,2);
%         ANorm = zeros(this.na,this.nb,this.nv/2,2);
%         A(:,:,:,1) = reshape(y1,[this.na this.nb this.nv/2]);
%         A(:,:,:,2) = reshape(y2,[this.na this.nb this.nv/2]);
%         ANorm(:,:,:,1) = reshape(y1Norm,[this.na this.nb this.nv/2]);
%         ANorm(:,:,:,2) = reshape(y2Norm,[this.na this.nb this.nv/2]);        
        y = DataTypes.CTData(this,A,ANorm,object.L);
    end
end