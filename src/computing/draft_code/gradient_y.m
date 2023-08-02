% Gradient discret vertical
function [nabla_y] = gradient_y(u)
[M,N]=size(u);
nabla_y(1:M-1,:)=u(2:M,:)-u(1:M-1,:);
nabla_y(M,:)=zeros(1,N);
end