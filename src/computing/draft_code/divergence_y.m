% Opérateur de divergence discrète verticale
function [div_y] = divergence_y(p2)
[M,N]=size(p2);

div_y(1,:)=p2(1,:);
div_y(2:M-1,:)=p2(2:M-1,:)-p2(1:M-2,:);
div_y(M,:)=-p2(M-1,:);
end
