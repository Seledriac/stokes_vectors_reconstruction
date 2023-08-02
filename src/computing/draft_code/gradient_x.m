% Gradient discret horizontal
function [nabla_x] = gradient_x(u)
[M,N]=size(u);
nabla_x(:,1:N-1)=u(:,2:N)-u(:,1:N-1);
nabla_x(:,N)=zeros(M,1);
end

