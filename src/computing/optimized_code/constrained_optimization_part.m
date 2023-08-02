% Renvoie le point proximal de t par rapport
% au terme d'attache aux données + l'indicatrice de la condition d'admissibilité en 1 pixel
% Résultat obtenu par conditions de KKT sur le problème de minimisation
% de l'opérateur proximal
function sol_ex = constrained_optimization_part(b,t,tau,mu)

alpha=mu+(1/tau);
beta=(mu/2)+(1/tau);
% variables intermédiaires
z0=(1/alpha)*(mu*b(1)+(1/tau)*t(1));
z1=(1/beta)*(mu*b(2)+(1/tau)*t(2));
z2=(1/beta)*(mu*b(3)+(1/tau)*t(3));

if(z0>=0)
    
    if(sqrt(z1*z1+z2*z2)<=z0)
        s0_ast=z0;
        s1_ast=z1;
        s2_ast=z2;        
    else
        s0_ast=(alpha/(alpha+beta))*z0+(beta/(alpha+beta))*sqrt(z1*z1+z2*z2);
        s1_ast=((beta/(alpha+beta))+(alpha/(alpha+beta))*(z0/sqrt(z1*z1+z2*z2)))*z1;
        s2_ast=((beta/(alpha+beta))+(alpha/(alpha+beta))*(z0/sqrt(z1*z1+z2*z2)))*z2;        
    end
    
else
    
    if(beta*sqrt(z1*z1+z2*z2)<=-alpha*z0)       
        s0_ast=0;
        s1_ast=0;
        s2_ast=0;        
    else    
        s0_ast=(alpha/(alpha+beta))*z0+(beta/(alpha+beta))*sqrt(z1*z1+z2*z2);
        s1_ast=((beta/(alpha+beta))+(alpha/(alpha+beta))*(z0/sqrt(z1*z1+z2*z2)))*z1;
        s2_ast=((beta/(alpha+beta))+(alpha/(alpha+beta))*(z0/sqrt(z1*z1+z2*z2)))*z2;
    end

end

sol_ex=[s0_ast;s1_ast;s2_ast];
