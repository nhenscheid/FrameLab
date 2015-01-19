function y=apply(this,object)
    this.setPara(object); %set the Gao parameter struct
    X0 = object.dataArray(:);
    if gpuDeviceCount==0
        if(this.verbose==true)
            disp('Computing forward cone beam transform with CPU');
            this.para
        end
        this.para.GPU = uint32(0);
        this.checkInputs(X0);
        y = Ax_fan_mf_cpu(X0,this.para);
    elseif gpuDeviceCount>0
        if(this.verbose==true)
            disp('Computing forward cone beam transform with GPU');
            this.para
        end
        this.para.GPU = uint32(1);
        this.checkInputs(X0);
        
        y = Ax_fan_mf(X0,this.para);
    end
    
    
    y = DataTypes.CTData(this,reshape(y,[this.nd,this.nv]),[],object.L);
end