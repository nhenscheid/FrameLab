function y=doScan(this,object)
    this.setPara(object); %set the Gao parameter struct
    X0 = object.dataArray(:);
    this.checkInputs(X0);
    size(X0)
    this.para
    if gpuDeviceCount==0
        disp('Computing forward cone beam transform with CPU');
        y = Ax_cone_mf(X0,this.para);
    elseif gpuDeviceCount>0
        disp('Computing forward cone beam transform with GPU');
        y = Ax_cone_mf_gpu(X0,this.para);
    end
end