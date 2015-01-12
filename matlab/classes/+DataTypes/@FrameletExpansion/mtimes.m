function c = mtimes(a,b)
    % multiply the framelet array b by the scalar a
    % 
    
    dim = b.dim;
    level = b.nLevel;
    type = b.frameletType; %Assuming both are of the same type?
    c = DataTypes.FrameletExpansion(dim,type,level,b.frameletArray);
    
    switch dim
        case 2
            [nx,ny] = size(b.frameletArray{1});
            for ki=1:level
                for ji = 1:nx
                    for jj=1:ny
                        c.frameletArray{ki}{ji,jj} = a*b.frameletArray{ki}{ji,jj};
                    end
                end
            end
           
            
        case 3
            [nx,ny,nz] = size(b.frameletArray{1});
            for ki=1:level
                for ji=1:nx
                    for jj=1:ny
                        for jk=1:nz
                            c.frameletArray{ki}{ji,jj,jk} = a*b.frameletArray{ki}{ji,jj,jk};
                        end
                    end
                end
            end
    end
    

end