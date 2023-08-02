
% Récupération des images de radiance
fI0=matfile(strcat(data_source,'I0.mat'));
I0=fI0.I0;
fI90=matfile(strcat(data_source,'I90.mat'));
I90=fI90.I90;
fI45=matfile(strcat(data_source,'I45.mat'));
I45=fI45.I45;
fI135=matfile(strcat(data_source,'I135.mat'));
I135=fI135.I135;


% Algorithme CP sans refitting pour reconstruction des vecteurs de Stokes
% à partir de l'image polarimétrique de synthèse bruitée
% Channels de variation totale coupl�s

M=size(I0,1);
N=size(I0,2);
MN = M*N;

% variables primales
S_hat=zeros(M,N,3);
S_new_hat=zeros(M,N,3);
S_bar_hat=S_hat;

% variables duales
p_hat=zeros(M,N,6);

% paramètres A et b pour le prox de S_hat
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
    
    % Ex�cution de l'it�ration
    [p_hat,S_new_hat]=cm_iteration_bi(norm1,norm2,nsup11,nsup12,nsup1r1,nsup1r2,S_hat,z_hat,M,N,MN,tau_bi,mu_bi,gamma_bi,b);
    
    % Accélération de la convergence en continuant dans la direction
    % donnée par la mise à jour de S_hat
    S_bar_hat=S_new_hat+theta_bi*(S_new_hat-S_hat);    
   
    S_hat=S_new_hat;
    
end

% Mesure du temps de calcul
toc


% Sauvegarde des variables
save(strcat(data_dest,'S_hat_mu_bi_',sprintf('%.1f',mu_bi),'.mat'),'S_hat');
