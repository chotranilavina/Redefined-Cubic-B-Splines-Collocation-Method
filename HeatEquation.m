clc
clear
%Equation: du/dt=C*d2u/dx2
%c=1
h = 0.05; x0 = 0; x1 = 1; b = 1; t0 = 0; t1 = 0.1; k = 0.01;
x = x0: h: x1;
sizex = x;
t = t0: k: t1;
M = length(x);
N = length(t);
p0 = sin(2*pi*x);
u = zeros(N,M);
u(1,:) = p0;
% M = M - 2;
for i = 2: N
[p_i] = fun_bspline(p0,h,k,M);
u(i,:) = p_i;
p0 = p_i;
end
V = zeros(N,M);
V(1,:) = sin(2*pi*x);
for i = 1: N
for j = 1: M
V(i,j) = exp(-4*t(i)*pi^2)*sin(2*pi*x(j));
U_exact(j) = exp((-4)*t(i)*pi^2)*sin(2*pi*x(j));
end
end
error = max(abs(u-V));
yerr = abs(u-V);
y=[x' p_i' U_exact' error']
%-----------------
figure(1)
surf(x,t,V)
%-----------------
% figure(1)
% plot(x,yerr)
%-----------------
figure(2)
surf(x,t,u)
%-----------------
% surf(x,t,u)
function YYp = fun_bspline( p0,h,k,M)
gammaj = 12*k + 4*h^2;
alphaj = -6*k + h^2;
betaj = -6*k + h^2;
A = zeros(M,M);
Mvalue = M;
A = diag([36*k ones(1,M-2)*gammaj 36*k],0) + diag([0 ones(1,M-2)*alphaj],1) + diag([ones(1,M-2)*betaj 0],-1);
ft = p0;
alpha0 = 0;
alpha1 = 0;
yp0 = alpha0;
yp1 = alpha1;
dm = h^2*ft;
dm(1) = dm(1) - yp0*(-6*k+h^2);
dm(M) = dm(M) - yp1*(-6*k+h^2);
Amatrix = A;
dm = dm';
sizeA = size(A);
sizedm = size(dm);
cc = A\dm;
ccm1 = 0 - 4*cc(1) - cc(2);
ccmp1 = 0 - 4*cc(M) - cc(M-1);
ccc = [ccm1 cc' ccmp1];
YYp = ccc(1:M) + 4*ccc(2:M+1) + ccc(3:M+2);
end
