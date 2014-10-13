function y=apply(this,object)
    this.setPara(object); %set the Gao parameter struct
    X0 = object.dataArray(:);
    size(X0)
    if gpuDeviceCount==0
        disp('Computing forward cone beam transform with CPU');
        this.para.GPU = uint32(0);
        this.checkInputs(X0);
        this.para
        y = Ax_fan_mf(X0,this.para);
    elseif gpuDeviceCount>0
        disp('Computing forward cone beam transform with GPU');
        this.para.GPU = uint32(1);
        this.checkInputs(X0);
        this.para
        y = Ax_fan_mf(X0,this.para);
    end
end