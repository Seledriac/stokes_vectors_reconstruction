
% Récupération des images de radiance
fI0=matfile(strcat(data_source,'I0.mat'));
I0=fI0.I0;
fI90=matfile(strcat(data_source,'I90.mat'));
I90=fI90.I90;
fI45=matfile(strcat(data_source,'I45.mat'));
I45=fI45.I45;
fI135=matfile(strcat(data_source,'I135.mat'));
I135=fI135.I135;


% Algorithme CP avec refitting pour reconstruction des vecteurs de Stokes
% à partir de l'image polarimétrique de synthèse bruitée

M=size(I0,1);
N=size(I0,2);
MN = M*N;

% variables primales, problème à solution biaisée
S_hat=zeros(M,N,3);
S_new_hat=zeros(M,N,3);
S_bar_hat=S_hat;

% variables primales, problème de refitting
S_til=zeros(M,N,3);
S_new_til=zeros(M,N,3);
S_bar_til=S_til;

% variables duales, problème à solution biaisée
p_hat=zeros(M,N,6);

% variables duales, problème de refitting
p_til=zeros(M,N,6);

% paramètres A et b pour le prox de S_hat et S_til
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
b=A'*[I0(:)';I90(:)';I45(:)';I135(:)'];


% Mesure du temps de calcul
tic

for k=1:Nb_iter_max

    % variables intermédiaires z_hat, argument du premier prox et utile
    % pour le co-support
    z_hat=p_hat+sigma_bi*gamma_bi*gradient(S_bar_hat);

    % récupération du co-support
    norm1=vecnorm(z_hat(:,:,1:2),2,3); % matrice des normes R2 sur le canal z_0 pixel-wise
    nsup11=norm1>1; % indices où norme > 1 matrice 2D
    nsup1r1=repmat(nsup11,1,1,2); % indices où norme > 1 matrice 3D
    norm2=vecnorm(z_hat(:,:,3:6),2,3); % matrice des normes R4 sur les canaux (z_1,z_2) pixel-wise
    nsup12=norm2>1; % indices où norme > 1 matrice 2D
    nsup1r2=repmat(nsup12,1,1,6); % indices où norme > 1 matrice 3D
    nsup1r2(:,:,1:2)=0;
    
    % Exécution des deux itérations en parallèle
%     f1=parfeval(pool,@cm_iteration_bi,2,norm1,norm2,nsup11,nsup12,nsup1r1,nsup1r2,S_hat,z_hat,p_hat,M,N,MN,tau_bi,mu_bi,gamma_bi,b);
%     f2=parfeval(pool,@cm_iteration_re,2,p_til,S_til,S_bar_til,sigma_re,nsup1r1,nsup1r2,z_hat,sigma_bi,tau_re,gamma_re,mu_re,M,N,MN,b);
%     [p_hat,S_new_hat]=fetchOutputs(f1);
%     [p_til,S_new_til]=fetchOutputs(f2);
    [p_hat,S_new_hat]=cm_iteration_bi(norm1,norm2,nsup11,nsup12,nsup1r1,nsup1r2,S_hat,z_hat,M,N,MN,tau_bi,mu_bi,gamma_bi,b);
    [p_til,S_new_til]=cm_iteration_re(p_til,S_til,S_bar_til,sigma_re,nsup1r1,nsup1r2,z_hat,sigma_bi,tau_re,gamma_re,mu_re,M,N,MN,b);
    
    % Accélération de la convergence pour le problème à solution biaisée en continuant dans la direction
    % donnée par la mise à jour de S_hat
    S_bar_hat=S_new_hat+theta_bi*(S_new_hat-S_hat);
    
    % Accélération de la convergence pour le problème de refitting en continuant dans la direction
    % donnée par la mise à jour de S_til  
    S_bar_til=S_new_til+theta_re*(S_new_til-S_til);
   
    S_hat=S_new_hat;
    S_til=S_new_til;
    
end

% Mesure du temps de calcul
toc


% Sauvegarde des variables
save(strcat(data_dest,'S_hat_mu_bi_',sprintf('%.1f',mu_bi),'.mat'),'S_hat');
save(strcat(data_dest,'S_til_mu_bi_',sprintf('%.1f',mu_bi),'_mu_re_',sprintf('%.1f',mu_re),'.mat'),'S_til');
