function checkInputs(this,X0)
    validateattributes(this.para.SO,{'single'},{'positive'});
    validateattributes(this.para.OD,{'single'},{'positive'});
    validateattributes(this.para.scale,{'single'},{'positive'});
    validateattributes(this.para.nx,{'uint32'},{'positive'});
    validateattributes(this.para.ny,{'uint32'},{'positive'});
    validateattributes(this.para.nz,{'uint32'},{'positive'});
    validateattributes(this.para.nt,{'uint32'},{'positive'});%\geq 1
    %validateattributes(this.para.sd_phi,{'single'},{'numel',this.nv});
    %validateattributes(this.para.sd_z,{'single'},{'numel',this.nv});
    validateattributes(this.para.y_det,{'single'},{'numel',this.na});
    validateattributes(this.para.z_det,{'single'},{'numel',this.nb});
    validateattributes(this.para.cos_phi,{'single'},{});
    validateattributes(this.para.sin_phi,{'single'},{});
    validateattributes(this.para.dy_det,{'single'},{'positive'});
    validateattributes(this.para.y_os,{'single'},{});
    validateattributes(this.para.dz_det,{'single'},{'positive'});
    validateattributes(this.para.dz,{'single'},{});
    validateattributes(this.para.id_X,{'uint32'},{});
    validateattributes(this.para.id_Y,{'uint32'},{});
    validateattributes(this.para.Nv,{'uint32'},{});
    validateattributes(this.para.tmp_size,{'uint32'},{});
    validateattributes(this.para.nv_block,{'uint32'},{});
    
    validateattributes(X0,{'single'},{});
end