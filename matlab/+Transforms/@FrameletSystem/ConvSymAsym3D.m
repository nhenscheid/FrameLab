function C=ConvSymAsym3D(A,M,b,L)

% function C=ConvSymAsym(A,M,b,L)
% The function implements 1D convolution with symmetric or antisymmetric
% boundary condition. The data are extened row by row, not block by block.
% For example, for the vector [1 2 3 4 5 6 7 8] and L=2, b='s', it is extended
% like [...4 3 2 1|1 2 3 4 5 6 7 8|8 7 6 5...], 
% not like [ 3 4:1 2|1 2 3 4 5 6 7 8|7 8:5 6....]
% A is the data, M is the filter, b is the boundary condition, L is the level.

% Written by Jian-Feng Cai. Modified to 3D by Bin Dong

[m,n]=size(A);
nM=length(M);
%n2M=floor(nM/2);
step=2^(L-1);

%====another implementation=========================
ker=zeros(step*(nM-1)+1,1);
ker(1:step:step*(nM-1)+1,1)=M;
lker=floor(length(ker)/2);

if b=='c'
    C=imfilter(A,ker,'circular');
else
    if b=='s'
        Ae=padarray(A,lker,'symmetric','both');
    elseif b=='a'
        Ae=padarray(A,lker,'symmetric','both');
        Ae(1:lker,:)=-Ae(1:lker,:);
        Ae(m+lker+1:m+2*lker,:)=-Ae(m+lker+1:m+2*lker,:);
    end
    C=conv2(Ae,ker,'valid');
end