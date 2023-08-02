% Opérateur de divergence discrète horizontale
function [div_x] = divergence_x(p1)
[M,N]=size(p1);

div_x(:,1)=p1(:,1);
div_x(:,2:N-1)=p1(:,2:N-1)-p1(:,1:N-2);
div_x(:,N)=-p1(:,N-1);
end

