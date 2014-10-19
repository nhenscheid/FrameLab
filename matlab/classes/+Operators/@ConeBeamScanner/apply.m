function y=apply(this,object)
    this.setPara(object); %set the Gao parameter struct
    X0 = object.dataArray(:);
    if gpuDeviceCount==0
        this.para.GPU = uint32(0);
        this.checkInputs(X0);
        if(this.verbose)
            disp('Computing forward cone beam transform with CPU');
            this.para
        end
        y = Ax_cone_mf(X0,this.para);
    elseif gpuDeviceCount>0
        this.para.GPU = uint32(1);
        this.checkInputs(X0);
        if(this.verbose)
            disp('Computing forward cone beam transform with GPU');
            this.para
        end
        y = Ax_cone_mf(X0,this.para);
    end
    
    y = DataTypes.ConeBeamData(reshape(y,[this.na this.nb this.nv]),this.para,object.L,'default');
end