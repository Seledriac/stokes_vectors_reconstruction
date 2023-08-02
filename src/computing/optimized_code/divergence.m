% Opérateur de divergence discrète de Y^k dans X^k
function [div] = divergence(p)
[M,N,k2]=size(p);
k = k2/2;
div=zeros(M,N,k);
div(:,2:N-1,:)=div(:,2:N-1,:)+p(:,2:N-1,1:2:k2)-p(:,1:N-2,1:2:k2);
div(:,1,:)=div(:,1,:)+p(:,1,1:2:k2);
div(:,N,:)=div(:,N,:)-p(:,N-1,1:2:k2);
div(2:M-1,:,:)=div(2:M-1,:,:)+p(2:M-1,:,2:2:k2)-p(1:M-2,:,2:2:k2);
div(1,:,:)=div(1,:,:)+p(1,:,2:2:k2);
div(M,:,:)=div(M,:,:)-p(M-1,:,2:2:k2);
end