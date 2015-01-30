% Kron test, 4D
clear all; clc;
f = @(a,b,t,z)sin(a).*cos(b).*exp(t).*z;
fa = @(a,b,t,z)cos(a).*cos(b).*exp(t).*z;
fb = @(a,b,t,z)-sin(a).*sin(b).*exp(t).*z;
faz = @(a,b,t,z)cos(a).*cos(b).*exp(t);

na = 128;
nb = 128;
nt = 128;
nz = 3; 

da = 1/(na-1);
db = 1/(nb-1);
dt = 1/(nt-1);
dz = 1/(nz-1);

a = 0:da:1;
b = 0:db:1;
t = 0:dt:1;
z = 0:dz:1;

[A,B,T,Z] = ndgrid(a,b,t,z);
u = f(A,B,T,Z);
ua_e = fa(A,B,T,Z);
ub_e = fb(A,B,T,Z);
uaz_e = faz(A,B,T,Z);

disp('Testing a derivative')
ea = ones(na,1);
za = zeros(na,1); 
Da = spdiags([-ea za ea],0:2,na-2,na);
tic
D2 = kron(eye(nb),Da);
D3 = kron(eye(nt),D2);
D = kron(eye(nz),D3)/(2*da);
toc

tic
ua = reshape(D*u(:),[na-2,nb,nt,nz]);
toc

ua_e2 = ua_e(2:end-1,:,:,:);
norm(ua(:)-ua_e2(:))/norm(ua_e2(:))

disp('Testing b derivative')
eb = ones(nb,1);
zb = zeros(nb,1); 
Db = spdiags([-eb zb eb],0:2,nb-2,nb);
D2 = kron(Db,eye(na));
D3 = kron(eye(nt),D2);
D = kron(eye(nz),D3)/(2*db);

ub = reshape(D*u(:),[na,nb-2,nt,nz]);

ub_e2 = ub_e(:,2:end-1,:,:);
norm(ub(:)-ub_e2(:))/norm(ub_e2(:))

disp('Testing az derivative')
ez = ones(nz,1);
zz = zeros(nz,1);
Da = spdiags([-ea za ea],0:2,na-2,na)/(2*da);
Dz = spdiags([-ez ez],0:1,nz-1,nz)/dz;
 
D1 = kron(eye(nb),Da);
D2 = kron(Dz,eye(nt));
D = kron(D2,D1);

uaz = reshape(D*u(:),[na-2,nb,nt,nz-1]);
uaz_e2 = uaz_e(2:end-1,:,:,1:end-1);
norm(uaz(:)-uaz_e2(:))/norm(uaz_e2(:))


