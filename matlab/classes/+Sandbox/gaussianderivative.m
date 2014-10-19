N=10000;
x = linspace(-pi,pi,N);
y = sin(x);
plot(x,y)

sig = 0.2;
G = exp(-x.^2/sig^2)/(sig*sqrt(pi));

figure; plot(x,G);

z = conv(y,G,'same')*2*pi/(N-1);

figure;plot(z)

norm(abs(z-y))/norm(y)


Gp = exp(-x.^2/sig^2).*(-2*x/sig^2)/(sig*sqrt(pi));
plot(x,Gp)

zp = conv(y,Gp)*2*pi/(N-1);
zpe = cos(x);

figure;plot(zp);
norm(abs(zp-zpe))/norm(zpe)