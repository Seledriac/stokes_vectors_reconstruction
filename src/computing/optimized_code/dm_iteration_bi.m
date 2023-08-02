function [p_new_hat,S_new_hat]=dm_iteration_bi(norm1,norm2,norm3,nsup11,nsup12,nsup13,nsup1r1,nsup1r2,nsup1r3,S_hat,z_hat,M,N,MN,tau_bi,mu_bi,gamma_bi,b)

    p_new_hat=z_hat;
    S_new_hat=zeros(M,N,3);

    % ########## UPDATE p_hat ##########

    % projection sur la boule de rayon 1 dans R^2 pour z_0_hat
    p_new_hat(nsup1r1)=p_new_hat(nsup1r1)./repmat(norm1(nsup11),2,1,1);
    % projection sur la boule de rayon 1 dans R^2 pour z_1_hat
    p_new_hat(nsup1r2)=p_new_hat(nsup1r2)./repmat(norm2(nsup12),2,1,1);
    % projection sur la boule de rayon 1 dans R^2 pour z_2_hat
    p_new_hat(nsup1r3)=p_new_hat(nsup1r3)./repmat(norm3(nsup13),2,1,1);
    
    % ########## UPDATE S_hat ##########

    % argument du 2ème prox    
    s_arg_hat=S_hat+tau_bi*gamma_bi*(divergence(p_new_hat));    
    s=[reshape(s_arg_hat(:,:,1),[1,MN]);reshape(s_arg_hat(:,:,2),[1,MN]);reshape(s_arg_hat(:,:,3),[1,MN])];

    % résultat du prox    
    S0_new_hat=zeros([1 MN]);
    S1_new_hat=zeros([1 MN]);
    S2_new_hat=zeros([1 MN]);
    
    % calcul du prox séparé : pixel par pixel
    for i=1:MN
        % Calcul du deuxième prox au pixel i :
        % fonction attache aux données + indicatrice de la contrainte
        % d'admissibilité physique
        sol_ex=constrained_optimization_part(b(:,i),s(:,i),tau_bi,mu_bi);        
        S0_new_hat(i)=sol_ex(1);
        S1_new_hat(i)=sol_ex(2);
        S2_new_hat(i)=sol_ex(3);
    end    
    S_new_hat(:,:,1)=reshape(S0_new_hat,[M,N]);
    S_new_hat(:,:,2)=reshape(S1_new_hat,[M,N]);
    S_new_hat(:,:,3)=reshape(S2_new_hat,[M,N]);

end
