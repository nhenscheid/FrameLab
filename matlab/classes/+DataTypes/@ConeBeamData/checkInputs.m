function checkInputs(this,dataArray,para,scanType)
%Check construction inputs for ConeBeamData 

    validateattributes(dataArray,{'single'},{'size',[length(para.y_det),length(para.z_det),para.Nv]});


end

