function c = plus(a,b)
    % Addition of framelet arrays.  a and b are of type FrameletExpansion.
    % 
    if(a.dim~=b.dim)
        error(['warning: inconsistent dimensions: a.dim = ',num2str(a.dim),'b.dim = ',num2str(b.dim)]);
    end
    
    dim = a.dim;
    level = a.nLevel;
    type = a.frameletType; %Assuming both are of the same type?
    c = DataTypes.FrameletExpansion(dim,type,level,a.frameletArray);
    
    switch dim
        case 2
            [nx,ny] = size(a.frameletArray{1});
            for ki=1:level
                for ji = 1:nx
                    for jj=1:ny
                        c.frameletArray{ki}{ji,jj} = a.frameletArray{ki}{ji,jj}+b.frameletArray{ki}{ji,jj};
                    end
                end
            end
           
            
        case 3
            [nx,ny,nz] = size(a.frameletArray{1});
            for ki=1:level
                for ji=1:nx
                    for jj=1:ny
                        for jk=1:nz
                            c.frameletArray{ki}{ji,jj,jk} = a.frameletArray{ki}{ji,jj,jk}+b.frameletArray{ki}{ji,jj,jk};
                        end
                    end
                end
            end
    end
    

end