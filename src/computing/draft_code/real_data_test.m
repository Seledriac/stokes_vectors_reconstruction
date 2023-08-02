%clear all
%close all
clc

% Image polarimétrique à 4 intensités par "superpixel"
I=imread('testfile_polarimetrie.tiff');

% Normalisation des intensités
I=double(I);
I=I/max(max(I));

% Récupération des intensités présentes dans les superpixels

I45=I(1:2:end,1:2:end);
%figure;imshow(I45,[]);
%size(I45)

I135=I(2:2:end,2:2:end);
%figure;imshow(I135,[]);
%size(I135)

I0=I(1:2:end,2:2:end);
%figure;imshow(I0,[]);
%size(I0)

I90=I(2:2:end,1:2:end);
%figure;imshow(I90,[]);
%size(I90)

% Les superpixels étant carrés, les dimensions des images
% de chaque canal sont identiques
[M,N]=size(I0);

% Action du filtre polarisant orienté à 4 angles sur Sin
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Moindres carrés classiques %%%%%%%%%%%%%%%%%%

for i=1:M
    for j=1:N

        % A^t * I
        b=A'*[I0(i,j);I90(i,j);I45(i,j);I135(i,j)];

        % pseudo-inverse = inv(A^t*A)
        sol_ex=inv(A'*A)*b;

        S0_least_squares(i,j)=sol_ex(1);
        S1_least_squares(i,j)=sol_ex(2);
        S2_least_squares(i,j)=sol_ex(3);
    
    end
end
clear sol_ex;

% 
% figure;imshow(S0_least_squares,[]);title('Composante S0 par les moindres carrés usuels');
% pause(2)
% figure;imshow(S1_least_squares,[]);title('Composante S1 par les moindres carrés usuels');
% pause(2)
% figure;imshow(S2_least_squares,[]);title('Composante S2 par les moindres carrés usuels');
% pause(2)
% figure;imshow(atan(S2_least_squares./S1_least_squares),[]);title('AOP par les moindres carrés usuels');
% pause(2);
% figure;imshow(sqrt(S1_least_squares.^2+S2_least_squares.^2)./S0_least_squares,[]);title('DOP par les moindres carré usuels');
% 
% [p,q]=find(S0_least_squares.^2<S1_least_squares.^2+S2_least_squares.^2)
% 
% for k=1:length(p)
%     
%    abs(S0_least_squares(p(k),q(k))^2- S1_least_squares(p(k),q(k))^2-S2_least_squares(p(k),q(k))^2) 
%     
%     
% end
% clear p q;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KKT Sans contrainte de régularité spatiale %%%%%%%%%%%%%%%%%%


% Résultat obtenu par conditions de KKT sur le problème de minimisation
% de l'attache au données + indicatrice de la condition d'admissibilité
% pour chaque pixel, le problème étant séparable par pixel.
for i=1:M
    for j=1:N
        
        b=A'*[I0(i,j);I90(i,j);I45(i,j);I135(i,j)];
        
        if(b(1)>=0)
            
            if(2*sqrt(b(2)*b(2)+b(3)*b(3))<=b(1))
                S0_no_regularization(i,j)=b(1);
                S1_no_regularization(i,j)=2*b(2);
                S2_no_regularization(i,j)=2*b(3);
                
            else
                S0_no_regularization(i,j)=(2/3)*(b(1)+sqrt(b(2)*b(2)+b(3)*b(3)));
                S1_no_regularization(i,j)=(2/3)*(1+b(1)/sqrt(b(2)*b(2)+b(3)*b(3)))*b(2);
                S2_no_regularization(i,j)=(2/3)*(1+b(1)/sqrt(b(2)*b(2)+b(3)*b(3)))*b(3);
                
            end
            
        else
            
            if(sqrt(b(2)*b(2)+b(3)*b(3))<=-b(1))
                S0_no_regularization(i,j)=0;
                S1_no_regularization(i,j)=0;
                S2_no_regularization(i,j)=0;
            else
                S0_no_regularization(i,j)=(2/3)*(b(1)+sqrt(b(2)*b(2)+b(3)*b(3)));
                S1_no_regularization(i,j)=(2/3)*(1+b(1)/sqrt(b(2)*b(2)+b(3)*b(3)))*b(2);
                S2_no_regularization(i,j)=(2/3)*(1+b(1)/sqrt(b(2)*b(2)+b(3)*b(3)))*b(3);
            end
                        
        end

    end
end

% [p,q]=find(S0_no_regularization.^2<S1_no_regularization.^2+S2_no_regularization.^2)
% for k=1:length(p)
%     
%    abs(S0_no_regularization(p(k),q(k))^2- S1_no_regularization(p(k),q(k))^2-S2_no_regularization(p(k),q(k))^2) 
%     
%     
% end
% clear p q;
% 
% figure;imshow(S0_no_regularization,[]);title('Composante S0 par la méthode KKT sans régularisation spatiale');
% pause(2)
% figure;imshow(S1_no_regularization,[]);title('Composante S1 par la méthode KKT sans régularisation spatiale');
% pause(2)
% figure;imshow(S2_no_regularization,[]);title('Composante S2 par la méthode KKT sans régularisation spatiale ');
% pause(2)
% figure;imshow(atan(S2_no_regularization./S1_no_regularization),[]);title('AOP par la méthode KKT sans régularisation spatiale');
% pause(2)
% figure;imshow(sqrt(S1_no_regularization.^2+S2_no_regularization.^2)./S0_no_regularization,[]);title('DOP par la méthode KKT sans régularisation spatiale');
% 


%%%%%%%%%%%%%%% KKT  avec régularisation spatiale couplée pour les canaux S1 et S2 + refitting %%%%%%%%%%%%%%%

% Dans le cas découplé, la contrainte se traduit par 
% \tau \sigma ||div ||^2<1 soit \tau \sigma <(1/8)
% avec la pondération par gamma, la condition de convergence s'exprime par
% tau sigma \gamma^2 <1/8

gamma_bi=1.0; % poids de TV dans la fonctionnelle du problème à solution biaisée
mu_bi=1.5; % poids du terme d'attache aux données dans la fonctionnelle du problème à solution biaisée

% sigma_bi * tau_bi * gamma_bi^2 = 1/16 < 1/8, convergence de CP pour problème à solution biaisée ok
sigma_bi=0.3;
tau_bi=1.0/(8.0*2*sigma_bi); 
theta_bi=1.0; % convergence ok

gamma_re=1.0; % poids de phi_SD dans la fonctionnelle du problème de refitting
mu_re=1.5; % poids du terme d'attache aux données dans la fonctionnelle du problème de refitting

% sigma_re * tau_re * gamma_re^2 = 1/16 < 1/8, convergence de CP pour problème de refitting ok
sigma_re=0.3;
tau_re=1.0/(8.0*2*sigma_re); 
theta_re=1.0; % convergence ok

M=size(I0,1);
N=size(I0,2);

% variables primales, problème à solution biaisée
S0_hat=zeros(M,N);
S1_hat=zeros(M,N);
S2_hat=zeros(M,N);
S0_new_hat=zeros(M,N);
S1_new_hat=zeros(M,N);
S2_new_hat=zeros(M,N);
S0_bar_hat=S0_hat;
S1_bar_hat=S1_hat;
S2_bar_hat=S2_hat;
S0_bar_new_hat=S0_hat;
S1_bar_new_hat=S1_hat;
S2_bar_new_hat=S2_hat;

% variables primales, problème de refitting
S0_til=zeros(M,N);
S1_til=zeros(M,N);
S2_til=zeros(M,N);
S0_new_til=zeros(M,N);
S1_new_til=zeros(M,N);
S2_new_til=zeros(M,N);
S0_bar_til=S0_til;
S1_bar_til=S1_til;
S2_bar_til=S2_til;

% variables duales, problème à solution biaisée
p1_S0_hat=zeros(M,N);
p2_S0_hat=zeros(M,N);
p1_S1_hat=zeros(M,N);
p2_S1_hat=zeros(M,N);
p1_S2_hat=zeros(M,N);
p2_S2_hat=zeros(M,N);

% variables duales, problème de refitting
p1_S0_til=zeros(M,N);
p2_S0_til=zeros(M,N);
p1_S1_til=zeros(M,N);
p2_S1_til=zeros(M,N);
p1_S2_til=zeros(M,N);
p2_S2_til=zeros(M,N);


Nb_iter_max=1000;
for k=1:Nb_iter_max


    % ################## ITÉRATION PROBLÈME A SOLUTION BIAISÉE ##################

    
    % ######### PREMIER PROX ##########

    % variables intermédiaires z_hat, argument du premier prox
    z1_S0_hat=p1_S0_hat+sigma_bi*gamma_bi*gradient_x(S0_bar_hat);
    z2_S0_hat=p2_S0_hat+sigma_bi*gamma_bi*gradient_y(S0_bar_hat);
    z1_S1_hat=p1_S1_hat+sigma_bi*gamma_bi*gradient_x(S1_bar_hat);
    z2_S1_hat=p2_S1_hat+sigma_bi*gamma_bi*gradient_y(S1_bar_hat);
    z1_S2_hat=p1_S2_hat+sigma_bi*gamma_bi*gradient_x(S2_bar_hat);
    z2_S2_hat=p2_S2_hat+sigma_bi*gamma_bi*gradient_y(S2_bar_hat);        

    % projection sur la boule de rayon 1 dans R^2 pour z_0_hat
    p1_S0_hat=z1_S0_hat./max(1,sqrt(z1_S0_hat.^2+z2_S0_hat.^2));
    p2_S0_hat=z2_S0_hat./max(1,sqrt(z1_S0_hat.^2+z2_S0_hat.^2)); 
    % projection sur la boule de rayon 1 dans R^4 pour z_tilde_hat
    p1_S1_hat=z1_S1_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));
    p2_S1_hat=z2_S1_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));   
    p1_S2_hat=z1_S2_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));
    p2_S2_hat=z2_S2_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));
   

    % ######### DEUXIÈME PROX ##########

    % argument du 2ème prox
    s0_arg_hat=S0_hat+tau*gamma_bi*(divergence_x(p1_S0_hat)+divergence_y(p2_S0_hat));
    s1_arg_hat=S1_hat+tau*gamma_bi*(divergence_x(p1_S1_hat)+divergence_y(p2_S1_hat));
    s2_arg_hat=S2_hat+tau*gamma_bi*(divergence_x(p1_S2_hat)+divergence_y(p2_S2_hat));
    
    % calcul du prox séparé : pixel par pixel
    for i=1:M
        for j=1:N

            % variable intermédiaire b pour le pixel (i,j)
            b=A'*[I0(i,j);I90(i,j);I45(i,j);I135(i,j)];
            % argument s du prox pour le pixel (i,j)
            s=[s0_arg_hat(i,j);s1_arg_hat(i,j);s2_arg_hat(i,j)];
            % Calcul du deuxième prox en (i,j) :
            % fonction attache aux données + indicatrice de la contrainte
            % d'admissibilité physique
            sol_ex=constrained_optimization_part(b,s,tau_bi,mu_bi);
            
            S0_new_hat(i,j)=sol_ex(1);
            S1_new_hat(i,j)=sol_ex(2);
            S2_new_hat(i,j)=sol_ex(3);

        end
    end
    
    % Accélération de la convergence pour le problème à solution biaisée en continuant dans la direction
    % donnée par la mise à jour de S_hat
    S0_bar_new_hat=S0_new_hat+theta_bi*(S0_new_hat-S0_hat);  
    S1_bar_new_hat=S1_new_hat+theta_bi*(S1_new_hat-S1_hat);  
    S2_bar_new_hat=S2_new_hat+theta_bi*(S2_new_hat-S2_hat);


    % ################## ITÉRATION PROBLÈME DE REFITTING ##################

    
    % ######### CALCUL DU CO-SUPPORT ##########
   
    % en chaque pixel (i,j), le canal 1 est pour la variable duale z_0
    % et le canal 2 pour les variables duales (z_1,Z_2)
    I_hat=zeros(M,N,2);
    for i=1:M
        for j=1:N
            % on remplace gamma_bi par 1 dans inégalité car
            % la projection se fait sur la boule unité dans
            % l'itération du problème à solution biaisée pour la variable duale
            if(sqrt(z1_S0_hat(i,j)^2+z2_S0_hat(i,j)^2) > 1)
                I_hat(i,j,1) = 1;
            end
            if(sqrt(z1_S1_hat(i,j)^2+z2_S1_hat(i,j)^2+z1_S2_hat(i,j)^2+z2_S2_hat(i,j)^2) > 1)
                I_hat(i,j,2) = 1;
            end
        end
    end


    % ######### PREMIER PROX ##########

    % variables intermédiaires z_til, argument du premier prox
    z1_S0_til=p1_S0_til+sigma_re*gradient_x(S0_bar_til);
    z2_S0_til=p2_S0_til+sigma_re*gradient_y(S0_bar_til);
    z1_S1_til=p1_S1_til+sigma_re*gradient_x(S1_bar_til);
    z2_S1_til=p2_S1_til+sigma_re*gradient_y(S1_bar_til);
    z1_S2_til=p1_S2_til+sigma_re*gradient_x(S2_bar_til);
    z2_S2_til=p2_S2_til+sigma_re*gradient_y(S2_bar_til);

    % calcul du prox séparé : pixel par pixel
    for i=1:M
        for j=1:N

            % Canal z_0
            z_til_1 = [
                z1_S0_til(i,j)
                z2_S0_til(i,j)
            ];
            if(I_hat(i,j,1) == 1) % pixel du co-support de grad(S_0)
                % variable intermédiaire estimant grad(S_0)
                nu_1 = [
                    z1_S0_hat(i,j)
                    z2_S0_hat(i,j)
                ];
                norm_nu_1 = sqrt(z1_S0(i,j)^2+z2_S0(i,j)^2);
                % on remplace encore gamma_bi par 1 car
                % la projection se fait sur la boule unité dans
                % l'itération du problème à solution biaisée pour la variable duale
                psi_1 = ( ( norm_nu_1 - 1 ) / ( sigma_bi * norm_nu_1 ) ) * nu_1;
                % calcul de prox_{phi_SD^*}(z_til_1,psi_1,I_hat)
                sol_ex_1 = refitting_part(z_til_1,psi_1);
            else % pixel hors du co-support                
                sol_ex_1 = z_til_1;
            end
            p1_S0_til(i,j)=sol_ex_1(1);
            p2_S0_til(i,j)=sol_ex_1(2);

            % Canaux (z_1,z_2)
            z_til_2 = [
                z1_S1_til(i,j)
                z2_S1_til(i,j)
                z1_S2_til(i,j)
                z2_S2_til(i,j)
            ];            
            if(I_hat(i,j,2) == 1) % pixel du co-support de (grad(S_1),grad(S_2))
                % variable intermédiaire estimant (grad(S_1),grad(S_2))
                nu_2 = [
                    z1_S1(i,j)
                    z2_S1(i,j)
                    z1_S2(i,j)
                    z2_S2(i,j)
                ];
                norm_nu_2 = sqrt(z1_S1(i,j)^2+z2_S1(i,j)^2+z1_S2(i,j)^2+z2_S2(i,j)^2);
                % on remplace encore gamma_bi par 1 car
                % la projection se fait sur la boule unité dans
                % l'itération du problème à solution biaisée pour la variable duale
                psi_2 = ( ( norm_nu_2 - 1 ) / ( sigma_bi * norm_nu_2 ) ) * nu_2;
                % calcul de prox_{phi_SD^*}(z_til_2,psi_2,I_hat)
                sol_ex_2 = refitting_part(z_til_2,psi_2,gamma_re);                    
            else % pixel hors du co-support
                sol_ex_2 = z_til_2;
            end
            p1_S1_til(i,j)=sol_ex_2(1);
            p2_S1_til(i,j)=sol_ex_2(2);
            p1_S2_til(i,j)=sol_ex_2(3);
            p2_S2_til(i,j)=sol_ex_2(4);
        end
    end


    % ######### DEUXIÈME PROX ##########

    % argument du 2ème prox
    s0_arg_til=S0_til+tau_re*(divergence_x(p1_S0_til)+divergence_y(p2_S0_til));
    s1_arg_til=S1_til+tau_re*(divergence_x(p1_S1_til)+divergence_y(p2_S1_til));
    s2_arg_til=S2_til+tau_re*(divergence_x(p1_S2_til)+divergence_y(p2_S2_til));

    % calcul du prox séparé : pixel par pixel
    for i=1:M
        for j=1:N            
            % variable intermédiaire b pour le pixel (i,j)
            b=A'*[I0(i,j);I90(i,j);I45(i,j);I135(i,j)];
            % argument du prox s pour le pixel (i,j)
            s=[s0_arg_til(i,j);s1_arg_til(i,j);s2_arg_til(i,j)];
            % Calcul du deuxième prox en (i,j) :
            % fonction attache aux données + indicatrice de la contrainte
            % d'admissibilité physique
            sol_ex=constrained_optimization_part(b,s,tau_re,mu_re);
            
            S0_new_til(i,j)=sol_ex(1);
            S1_new_til(i,j)=sol_ex(2);
            S2_new_til(i,j)=sol_ex(3);
        end
    end
    
    % Accélération de la convergence pour le problème de refitting en continuant dans la direction
    % donnée par la mise à jour de S_til
    S0_bar_til=S0_new_til+theta_bi*(S0_new_til-S0_til);  
    S1_bar_til=S1_new_til+theta_bi*(S1_new_til-S1_til);  
    S2_bar_til=S2_new_til+theta_bi*(S2_new_til-S2_til);
   
    S0_hat=S0_new_hat;
    S1_hat=S1_new_hat;
    S2_hat=S2_new_hat;
    S0_bar_hat=S0_bar_new_hat;
    S1_bar_hat=S1_bar_new_hat;
    S2_bar_hat=S2_bar_new_hat;
    S0_til=S0_new_til;
    S1_til=S1_new_til;
    S2_til=S2_new_til;
    
end

% Solution biaisée de l'algorithme CP appliqué à
% l'image polarimétrique d'entrée I
S0_KKT_coupled_regularization=S0_hat;
S1_KKT_coupled_regularization=S1_hat;
S2_KKT_coupled_regularization=S2_hat;

% Solution débiaisée de l'algorithme CP appliqué à
% l'image polarimétrique d'entrée I
S0_KKT_coupled_regularization_refitted=S0_til;
S1_KKT_coupled_regularization_refitted=S1_til;
S2_KKT_coupled_regularization_refitted=S2_til;

clear S0_hat S1_hat S2_hat S0_bar_hat S1_bar_hat S2_bar_hat S0_new_hat S1_new_hat S2_new_hat S0_bar_new_hat S1_bar_new_hat S2_bar_new_hat;
clear S0_til S1_til S2_til S0_bar_til S1_bar_til S2_bar_til S0_new_til S1_new_til S2_new_til S0_bar_new_til S1_bar_new_til S2_bar_new_til;

% Vérification de la condition d'admissibilité physique
[p,q]=find(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2>S0_KKT_coupled_regularization.^2)
for k=1:length(p)    
   abs(S0_KKT_coupled_regularization(p(k),q(k))^2- S1_KKT_coupled_regularization(p(k),q(k))^2-S2_KKT_coupled_regularization(p(k),q(k))^2)
end
clear p q;

% figure;imshow(atan(S2_KKT_coupled_regularization./S1_KKT_coupled_regularization),[]);
% figure;imshow(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization,[]);

% figure;subplot(221),imshow(S0_least_squares,[]);
% subplot(222),imshow(S0_no_regularization,[]);
% subplot(223),imshow(S0_KKT_coupled_regularization,[]);

% 
% 
% figure;subplot(221),imshow(S1_least_squares,[]);
% subplot(222),imshow(S1_no_regularization,[]);
% subplot(223),imshow(S1_KKT_coupled_regularization,[]);
% 
% 
% figure;subplot(221),imshow(S2_least_squares,[]);
% subplot(222),imshow(S2_no_regularization,[]);
% subplot(223),imshow(S2_KKT_coupled_regularization,[]);
% 
% 

% figure;subplot(221),imshow(atan(S2_least_squares./S1_least_squares),[]);
% subplot(222),imshow(atan(S2_no_regularization./S1_no_regularization),[]);
% subplot(223),imshow(atan(S2_KKT_coupled_regularization./S1_KKT_coupled_regularization),[]);

% 
% % figure;subplot(221),imshow(sqrt(S1_least_squares.^2+S2_least_squares.^2)./S0_least_squares,[]);
% % subplot(222),imshow(sqrt(S1_no_regularization.^2+S2_no_regularization.^2)./S0_no_regularization,[]);
% % subplot(223),imshow(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization,[]);
% % 
% figure;imshow(sqrt(S1_least_squares.^2+S2_least_squares.^2)./S0_least_squares,[]);
% figure;imshow(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization,[]);
% 
% figure;imshow(S0_least_squares,[]);
% figure;imshow(S0_KKT_coupled_regularization,[]);
% 
% figure;imshow(S1_least_squares,[]);
% figure;imshow(S1_KKT_coupled_regularization,[]);
% 
% figure;imshow(S2_least_squares,[]);
% figure;imshow(S2_KKT_coupled_regularization,[]);
% AOP_KKT_coupled_regularization=atan(S2_KKT_decoupled_regularization./S1_KKT_decoupled_regularization);