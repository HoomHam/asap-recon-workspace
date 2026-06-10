function [kernel,u] = createKBkernel(width,beta,length)
u = (0:length-1)/(length-1) * width/2;
f = beta*sqrt(1-(2*u/width).^2);
kernel = besseli(0,f)./width;
kernel = (kernel/max(kernel))';
u = u';
end