function S_new_hat=gma_iteration_bi(S_hat,M,N,MN,tau_bi,mu_bi,b)

    S_new_hat=zeros(M,N,3);

    % ########## UPDATE S_hat ##########

    % argument du prox    
    s_arg_hat=S_hat; % (car p^n+1 == 0)
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
