function [p_new_til,S_new_til] = gma_iteration_re(p_til,S_til,S_bar_til,sigma_re,nsup1r1,nsup1r2,z_hat,sigma_bi,tau_re,gamma_re,mu_re,M,N,MN,b)

    % ######### UPDATE P_TIL ##########

    % variables intermédiaires z_til, argument du premier prox
    z_til=p_til+sigma_re*gradient(S_bar_til);
    
    % reshape de z_til pour la boucle parfor
    z_til_1=[
        reshape(z_til(:,:,1),[1,MN]);
        reshape(z_til(:,:,2),[1,MN]);
        ];
    z_til_2=[
        reshape(z_til(:,:,3),[1,MN]);
        reshape(z_til(:,:,4),[1,MN]);
        reshape(z_til(:,:,5),[1,MN]);
        reshape(z_til(:,:,6),[1,MN]);
        ];    

    % reshape du co-support pour la boucle parfor
    nsup1r1_rshp=reshape(nsup1r1(:,:,1),[1,MN]);
    nsup1r2_rshp=reshape(nsup1r2(:,:,3),[1,MN]);
%     % duplication pour le calcul du prox sans boucle
%     nsup1r1_r=repmat(nsup1r1_rshp(1,:)==1,2,1);
%     nsup1r2_r=repmat(nsup1r2_rshp(1,:)==1,4,1);

    % variables intermédiaires psi
    nu_1 = [
        reshape(z_hat(:,:,1),[1,MN]);
        reshape(z_hat(:,:,2),[1,MN]);
        ];
    norm_nu_1 = repmat(vecnorm(nu_1,2,1),2,1); % norme par colonne
    % on remplace encore gamma_bi par 1 car
    % la projection se fait sur la boule unité dans
    % l'itération du problème à solution biaisée pour la variable duale
    psi_1=nu_1.*(norm_nu_1-1)./(sigma_bi*norm_nu_1);

    nu_2 = [
        reshape(z_hat(:,:,3),[1,MN]);
        reshape(z_hat(:,:,4),[1,MN]);
        reshape(z_hat(:,:,5),[1,MN]);
        reshape(z_hat(:,:,6),[1,MN]);
        ];
    norm_nu_2 = repmat(vecnorm(nu_2,2,1),4,1); % norme par colonne
    psi_2=nu_2.*(norm_nu_2-1)./(sigma_bi*norm_nu_2);

    % pour les pixels hors du co-support, initialisation à z_til
    p_new_til_1=z_til_1(1:2,:);
    p_new_til_2=z_til_2(1:4,:);

    % calcul du prox séparé : pixel par pixel    
    for i=1:MN
        if(nsup1r1_rshp(i)==1)
            p_new_til_1(:,i) = refitting_part(z_til_1(:,i),psi_1(:,i),gamma_re);
        end
        if(nsup1r2_rshp(i)==1)
            p_new_til_2(:,i) = refitting_part(z_til_2(:,i),psi_2(:,i),gamma_re);
        end
    end
%     % calcul du prox sans boucle
%     nb1=sum(nsup1r1_rshp,'all');
%     nb2=sum(nsup1r2_rshp,'all');    
%     if ~isempty(p_new_til_1(nsup1r1_r==1))
%         p_new_til_1(nsup1r1_r==1) = refitting_part_2( ...
%             reshape(z_til_1(nsup1r1_r==1),[2,nb1]), ...
%             reshape(psi_1(nsup1r1_r==1),[2,nb1]), ...
%             gamma_re,2);
%     end
%     if ~isempty(p_new_til_2(nsup1r2_r==1))
%         p_new_til_2(nsup1r2_r==1) = refitting_part_2( ...
%             reshape(z_til_2(nsup1r2_r==1),[4,nb2]), ...
%             reshape(psi_2(nsup1r2_r==1),[4,nb2]), ...
%             gamma_re,4);
%     end
    
    p_new_til=zeros(M,N,6);
    reshape(p_new_til_1(1,:),[M N]);
    p_new_til(:,:,1)=reshape(p_new_til_1(1,:),[M,N]);
    p_new_til(:,:,2)=reshape(p_new_til_1(2,:),[M,N]);
    p_new_til(:,:,3)=reshape(p_new_til_2(1,:),[M,N]);
    p_new_til(:,:,4)=reshape(p_new_til_2(2,:),[M,N]);
    p_new_til(:,:,5)=reshape(p_new_til_2(3,:),[M,N]);
    p_new_til(:,:,6)=reshape(p_new_til_2(4,:),[M,N]);

    % ######### UPDATE S_til ##########

    % argument du 2ème prox
    s_arg_til=S_til+tau_re*gamma_re*(divergence(p_new_til));
    s=[reshape(s_arg_til(:,:,1),[1,MN]);reshape(s_arg_til(:,:,2),[1,MN]);reshape(s_arg_til(:,:,3),[1,MN])];

    % résultat du prox
    S0_new_til=zeros([1 MN]);
    S1_new_til=zeros([1 MN]);
    S2_new_til=zeros([1 MN]);

    % calcul du prox séparé : pixel par pixel
    for i=1:MN
        % Calcul du deuxième prox au pixel i :
        % fonction attache aux données + indicatrice de la contrainte
        % d'admissibilité physique
        sol_ex=constrained_optimization_part(b(:,i),s(:,i),tau_re,mu_re);        
        S0_new_til(i)=sol_ex(1);
        S1_new_til(i)=sol_ex(2);
        S2_new_til(i)=sol_ex(3);
    end
    S_new_til=zeros(M,N,3);
    S_new_til(:,:,1)=reshape(S0_new_til,[M,N]);
    S_new_til(:,:,2)=reshape(S1_new_til,[M,N]);
    S_new_til(:,:,3)=reshape(S2_new_til,[M,N]);

end
