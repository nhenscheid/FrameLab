function u = ball3d(N)
    % Create a characteristic function of a unit ball in 3d
    dx = 2/N;
    x = dx*((-N/2:N/2-1)+1/2);
    y = x;
    z = x;
    u = zeros(N,N,N);

    disp('Creating 3d unit ball array');
    for i=1:N
        for j=1:N
            for k=1:N
                if( sqrt(x(i)^2+y(j)^2+z(k)^2) <=1)
                    u(i,j,k)=1;
                end
            end
        end
    end
    
end %ball3d
