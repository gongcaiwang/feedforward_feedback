r=[-10000:10000]/20000;
x=[-10000:10000]/20000;
sigma=0.5;
sigma2=0.1;
f=exp(-x.^2/sigma);
% plot(x,f);
xx =abs(x);
alpha = f;
% alpha=0;
% alpha = ones(size(alpha));


for ii=1:20001
    ff = ff + alpha(ii) * exp(-(x-r(ii)).^2/sigma).*exp(1i*phi(ii));
end

g = exp(-x.^2/sigma2);
loss = sum((g-abs(ff)).^2);

for 
grad = diff(loss)/diff(alpha);
alpha = alpha + 0.001*grad;

grad = diff(loss)/diff(phi);
phi = phi + 0.001*grad;
end