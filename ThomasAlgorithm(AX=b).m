%of the form AX=B
n=input('enter the order for the matrix');
for(i=1:n)
    for(j=1:n)
      a(i,j)=input('enter the element of coefficient matrix');
    end
end
for i=1:n
    r(i)=input('enter the RHS');
end
R(1)=0;
P=zeros(1,n);
Q=zeros(1,n-1);
R=zeros(1,n);
Y=zeros(1,n-1);
for i=1:n
      P(i)=a(i,i);
end
for i=1:n-1
    Q(i)=a(i,i+1);
end
for i=1:n-1
    R(i+1)=a(i+1,i);
end
Y(1)=Q(1)/P(1);
for i=2:n-1
    Y(i)=Q(i)/(P(i)-R(i)*Y(i-1));
end
W(1)=r(1)/P(1);
for i=2:n
    W(i)=(r(i)-R(i)*W(i-1))/(P(i)-R(i)*Y(i-1));
end
x(n)=W(n);
for i=n-1:-1:1
    x(i)=W(i)-Y(i)*x(i+1);
end