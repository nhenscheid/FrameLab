function checkInputs(this,dataArray,para,scanType)
%Check construction inputs for ConeBeamData 
    validatestring(scanType,{'cone','fan'});
    if(strcmp(scanType,'cone'))
        validateattributes(dataArray,{'single'},{'size',[length(para.y_det),length(para.z_det),para.Nv]});
    elseif(strcmp(scanType,'fan'))
        validateattributes(dataArray,{'single'},{'size',[length(para.y_det),para.Nv]});
    end
end

