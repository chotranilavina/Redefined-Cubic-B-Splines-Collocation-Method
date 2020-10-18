%ut=epsilon*uxx-u*ux
%u(x,0)=f(x)
%u(0,t)=u(1,t)=0
epsilon=1;
dt=0.1;
x=0:0.125:1;
N=length(x);
T=1;
M=T/dt;
n=16;
U=zeros(n,M);
a=-1*epsilon*dt;
U=sym('U', [n,M]);
m = n; 
dt = T/m;
alpha = 0;  beta = 0;  % BCs
sigma =  1/2;  
h = 1/n; 
i=2;
U=zeros(N,M);
K_exact=[0 0.135829 0.253638 0.336742 0.371577 0.350123 0.272582 0.149239 0];
for j=2:M
    b=dt*U(1,j);
    c=(1+dt*diff(U(1,j)));
    d=U(i,j-1)+dt*U(i-1,j)*diff(U(i-1,j));
    syms x Y(x);
    D1Y = diff(Y,x);
    D2Y = diff(Y,x,2);
    Eqn = a*D2Y + b*D1Y + c*Y == d;
    %yode = odeToVectorField(Eqn);
    %Yodefcn = @(x,Y) [y(2);-b.*Y(2)/a-c.*Y(1)/a+d/a];
    %Yodefcn = matlabFunction(yode, 'Vars',[x Y]);
    xspan = [0 1];
    DY=diff(Y);
    DY0 = [0 0];
    DY1 = [0 0];
    %[X,Y] = ode45(Yodefcn, xspan, Dy0,Dy1);
    %v=[y1;y2];
    %Dy=diff(y1,x);
    %S=dsolve(diff(v)==A*v+B, Dy(1)==0,Dy(0)==0);
    %pp = spline(X,Y(:,1));
    %disp(pp.coefs);
   % xcoord = X;
   % ycoord = Y;
   % order = 3;
   % p(:,1) = xcoord;
   % p(:,2) = ycoord(:,1);
   % n1 = size(p,1);
   % pspl = math(p,n1,order);
   % plot(pspl(:,1),pspl(:,2));
   % hold on;
   % plot(p(:,1),p(:,2));
   % U(i,j)=Y1(1);
end
err=0;
x=0:0.125:1;
N=length(x);
for t=0:10:100
    for k=1:N
            Z(k)=exp(-1*t*0.01)*sin(pi*x(k));
            if(t==10)
            err=max(err,Z(k)-K_exact(k));
            end
    end
end
err=(err/(9.78966*(n^2))*100*h);
% u_t-u_xx +(2+x+x^2)u_x+(1+x^2)*u(x,t)=(sin(pi*x*(1-x)))
% u(x,0) = 0; u(0,t) = 0; u(1,t) = 0
clear y;
clear x;
clear t;
m =n; 
n=N;
T = 1;
dt = T/m;
alpha = 0;  beta = 0;  % BCs
sigma =  1/2;  
h = 1/n; 
for i = 1:n+1
x(i) = (i-1)*h;
end
for j = 1:m+1
    t(j) = (j-1)*dt;
end

for j = 1:m+1
for i = 1:n+1
        a(i,j) = 2+x(i)+x(i).^2;
        b(i,j) = (1+x(i).^2)/2;
        f(i,j) = sin(pi*x(i).*(1-x(i)));
        r(i,j) = b(i,j)+2/dt;
        s(i,j) = b(i,j)-2/dt;
        E(i,j) = -6/h^2 - 3*a(i,j)/h + r(i,j);
        F(i,j) = 12/h^2 + 4*r(i,j);
        G(i,j) = -6/h^2 + 3*a(i,j)/h + r(i,j);
        H(i,j) = 6/h^2 + 3*a(i,j)/h - s(i,j);
        I(i,j) = -12/h^2 - 4*s(i,j);
        J(i,j) = 6/h^2 - 3*a(i,j)/h - s(i,j);
end
end
u = zeros(n+1,m+1);
for i = 1:n+1
    for j=1:m+1
    u(i, 1) = 0; 
    u(1,j) = 0;
    u(n+1,j) = 0;
end
end

%......................... Matrix A .......................................
j = 1;
for j = 1:m
A = zeros(n+1,n+1);
A(1,1) = F(1,j+1)-4*E(1,j+1); A(1,2) = G(1,j+1)-E(1,j+1);
for i  = 2:n
        A(i,i-1) = E(i,j+1); 
        A(i,i)   = F(i,j+1);
        A(i,i+1) = G(i,j+1);
end
A(n+1,n) = E(n+1,j+1)-G(n+1,j+1);  A(n+1,n+1) = F(n+1,j+1)-4*G(n+1,j+1);

%......................... Matrix B .......................................
B = zeros(n+1,n+1);
B(1,1) = I(1,j)-4*H(1,j); B(1,2) = J(1,j)-H(1,j);
for i  = 2:n
        B(i,i-1) = H(i,j); 
        B(i,i)   = I(i,j);
        B(i,i+1) = J(i,j);
end
B(n+1,n) = H(n+1,j)-J(n+1,j);  B(n+1,n+1) = I(n+1,j)-4*J(n+1,j);

%............. Matrix C .................................................
C = zeros(n+1,m+1);
C(1,j+1) = 2*f(1,j+1)-alpha*E(1,j+1);
for i = 2:n
     C(i,:) = 2*f(i,j+1);
end
C(n+1,j+1) = 2*f(n+1,j+1)-beta*G(n+1,j+1);
end
v = zeros(n+1,m+1);
for i = 2:n
for j = 2:m+1
    v(:,j) = inv(A)*B*v(:,j-1)+inv(A)*C(:,j);
    u(i,j) = v(i-1,j)+ 4*v(i,j)+v(i+1,j); 
 end
end

% % % % % ============ Double mesh ========
 m1 = 2*m;
 n1 = 2*n;
 dt1 = T/m1;
 
 for j = 1:m1+1
     t1(j) = (j-1)*dt1;
 end
 
x1(1) = 0; x1(2) = (x(1)+x(2))/2;
for i = 2:n
x1(2*i-1) = x(i);
x1(2*i) = (x(i+1)+x(i))/2;
end
x1(2*n+1) = x(n+1);

h1 = x1(n1)-x1(n1-1);
  
 
for j = 1:m1+1
for i = 1:n1+1
        a1(i,j) = 2+x1(i)+x1(i).^2;
        b1(i,j) = (1+x1(i).^2)/2;
        f1(i,j) = sin(pi*x1(i).*(1-x1(i)));
        r1(i,j) = b1(i,j)+2/dt1;
        s1(i,j) = b1(i,j)-2/dt1;
        E1(i,j) = -6/h1^2 - 3*a1(i,j)/h1 + r1(i,j);
        F1(i,j) = 12/h1^2 + 4*r1(i,j);
        G1(i,j) = -6/h1^2 + 3*a1(i,j)/h1 + r1(i,j);
        H1(i,j) = 6/h1^2 + 3*a1(i,j)/h1 - s1(i,j);
        I1(i,j) = -12/h1^2 - 4*s1(i,j);
        J1(i,j) = 6/h1^2 - 3*a1(i,j)/h1 - s1(i,j);
end
end
u1 = zeros(n1+1,m1+1);
for i = 1:n1+1
    for j=1:m1+1
    u1(i, 1) = 0; 
    u1(1,j) = 0;
    u1(n1+1,j) = 0;
end
end

%......................... Matrix A .......................................

j = 1;
for j = 1:m1
A1 = zeros(n1+1,n1+1);
A1(1,1) = F1(1,j+1)-4*E1(1,j+1); A1(1,2) = G1(1,j+1)-E1(1,j+1);
for i  = 2:n1
        A1(i,i-1) = E1(i,j+1); 
        A1(i,i)   = F1(i,j+1);
        A1(i,i+1) = G1(i,j+1);
end
A1(n1+1,n1) = E1(n1+1,j+1)-G1(n1+1,j+1);  A1(n1+1,n1+1) = F1(n1+1,j+1)-4*G1(n1+1,j+1);

% ......................... Matrix B .......................................
B1 = zeros(n1+1,n1+1);
B1(1,1) = I1(1,j)-4*H1(1,j); B1(1,2) = J1(1,j)-H1(1,j);
for i  = 2:n1
        B1(i,i-1) = H1(i,j); 
        B1(i,i)   = I1(i,j);
        B1(i,i+1) = J1(i,j);
end
B1(n1+1,n1) = H1(n1+1,j)-J1(n1+1,j);  B1(n1+1,n1+1) = I1(n1+1,j)-4*J1(n1+1,j);

%............. Matrix C .................................................
C1 = zeros(n1+1,m1+1);
C1(1,j+1) = 2*f1(1,j+1)-alpha*E1(1,j+1);
for i = 2:n1
     C1(i,:) = 2*f1(i,j+1);
end
C1(n1+1,j+1) = 2*f1(n1+1,j+1)-beta*G1(n1+1,j+1);
end
v1 = zeros(n1+1,m1+1);
for i = 2:n1
for j = 2:m1+1
    v1(:,j) = inv(A1)*B1*v1(:,j-1)+inv(A1)*C1(:,j);
    u1(i,j) = v1(i-1,j)+ 4*v1(i,j)+v1(i+1,j); 
 end
end

for i = 1:n+1
for j = 1:m+1
y(i,j) = u1(2*i-1, 2*j-1);                  
error(i,j) = abs(u(i,j)-y(i,j));   
end
end
%error=max(max(error))
mesh(x,t,u');