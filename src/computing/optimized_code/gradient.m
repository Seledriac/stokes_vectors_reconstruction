% Gradient discret de X^k dans Y^k
function [nabla] = gradient(u)
[M,N,k]=size(u);
nabla = zeros(M,N,2*k);
nabla(:,1:N-1,1:2:2*k)=u(:,2:N,:)-u(:,1:N-1,:);
nabla(:,N,1:2:2*k)=zeros(M,1,k);
nabla(1:M-1,:,2:2:2*k)=u(2:M,:,:)-u(1:M-1,:,:);
nabla(M,:,2:2:2*k)=zeros(1,N,k);
end