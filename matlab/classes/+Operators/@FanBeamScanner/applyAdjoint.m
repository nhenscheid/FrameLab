function Aty=applyAdjoint(this,object)

    if(isempty(fieldnames(this.para)))
        % If the scanner's parameters aren't set for some reason,
        % import para from the object
        this.para = object.para;
    else
        if(~structcmp(this.para,object.para))
            error('Object scan parameters and scanner parameters do not agree!')
        end
    end
    
    X0 = object.dataArray;
    if gpuDeviceCount==0
        if(this.verbose==true)
            disp('Computing adjoint fan beam transform with CPU');
            this.para
        end
        this.para.GPU = uint32(0);
        this.checkInputs(X0);
        y = Atx_fan_mf_cpu(X0,this.para);
    elseif gpuDeviceCount>0
        if(this.verbose==true)
            disp('Computing adjoint fan beam transform with GPU');
            this.para
        end
        this.para.GPU = uint32(1);
        this.checkInputs(X0);
        y = Atx_fan_mf(X0,this.para);
    end
    Aty = DataTypes.ObjectData(2,reshape(y,[this.para.nx this.para.ny]),object.L);
end