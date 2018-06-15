A = [3 -1 1; 0 2 1; -1 1 4];
b = [-1; 0; 1];
x0 = [1; 1; 0];
xopt = [-5/12; -1/12; 1/6];
[x,JError] = Jacobi (A, b, x0, xopt);
[x,GError]= GaussSeide(A, b, x0, xopt);

figure
hold on 
plot(1:20,JError,'b--o')
plot(1:20,GError,'r*')
legend('Gauss Seidel Method','Jacobi Method')
ylabel('Error')
xlabel('Number of iterations')
title('Jacobi Method and Gauss Seidel Method')
hold off

function [x,JError] = Jacobi (A, b, x0, xopt)
n=size(x0,1);
Jitr=0;
maxitr=20;
JError = [];
while Jitr < maxitr
    xold=x0;
    for i=1:n
        sigma=0;
        for j=1:n
            if j~=i
                sigma=sigma+A(i,j)*xold(j);
            end
        end
        x0(i)=(1/A(i,i))*(b(i)-sigma);
    end
    Jitr=Jitr+1;
    error=max(abs(xopt-x0));
    JError=[JError;error];
end
x = x0;
JError = JError;
end

function[x,GError]= GaussSeide(A, b, x0, xopt)
n=size(x0,1);
Gitr=0;
maxitr=20;
GError = [];
while Gitr < maxitr
    x_old=x0;
    for i=1:n
        sigma=0;
        for j=1:i-1
            sigma=sigma+A(i,j)*x0(j);
        end
        
        for j=i+1:n
            sigma=sigma+A(i,j)*x_old(j);
        end
        x0(i)=(1/A(i,i))*(b(i)-sigma);
    end
    Gitr=Gitr+1;
    error=max(abs(xopt-x0));
    GError=[GError;error];
end
x = x0;
GError = GError;
end
