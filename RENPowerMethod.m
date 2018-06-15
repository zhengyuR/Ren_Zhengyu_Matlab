A=hilb(3);
x0=[1;1;1];
nit=1000;
PowerMethod(A, x0, nit)
[V,D] = eig(A)

function PowerMethod(A, x0, nit)
x = x0;
N = size(A, 1);
V = zeros(N, N);
l = zeros(1, N);
for n = 1:N;
    for i = 1:nit;
    vnew= A*x;
    lambda = norm(vnew)/norm(x);
    x = vnew/norm(vnew); 
    end
    V(:,n)= x/norm(x);
    l(n) = lambda;
    A = A-lambda/norm(V(:,n))*V(:,n)*transpose(V(:,n));
end
V = V
D = diag(l)
condNumb = max(l)/min(l)
end