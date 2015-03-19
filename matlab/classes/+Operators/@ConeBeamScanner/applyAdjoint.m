function y=applyAdjoint(this,object,N)
    %***Apply adjoint cone beam operator to an object of type
    %    CTData('cone')
    %    N = [nx,ny,nz]  (in case we want the adjoint array to be a different
    %    size than the original array...)
    N = uint32(N);
    
    if(isempty(fieldnames(this.para)))
        this.para = object.para; %If the scanner's parameters aren't set,import para from the object
        
    else
        if(~structcmp(this.para,object.scanner.para))
            error('Object scan parameters and scanner parameters do not agree!')
        end
    end

    %X0 = object.dataArray(:);  I think the adjoint wants the full array?
    this.para.nx = N(1);
    this.para.ny = N(2);
    this.para.nz = N(3);
    X0 = object.dataArray;
    if (gpuDeviceCount==0 || this.GPU==0)
        this.para.GPU = uint32(0);
        this.checkInputs(X0);
        if(this.verbose)
            disp('Computing adjoint cone beam transform with CPU');
            this.para
        end
        y = Atx_cone_mf(X0,this.para);
    elseif (gpuDeviceCount>0 && this.GPU ==1)
        this.para.GPU = uint32(1);
        this.checkInputs(X0);
        if(this.verbose)
            disp('Computing adjoint cone beam transform with GPU');
            this.para
        end
        y = Atx_cone_mf(X0,this.para);
    end
    
    y = DataTypes.ObjectData(3,reshape(y,[this.para.nx this.para.ny this.para.nz]),object.L);
end