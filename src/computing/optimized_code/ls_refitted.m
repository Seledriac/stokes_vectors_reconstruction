
% Récupération des images de radiance
fI0=matfile(strcat(data_source,'I0.mat'));
I0=fI0.I0;
fI90=matfile(strcat(data_source,'I90.mat'));
I90=fI90.I90;
fI45=matfile(strcat(data_source,'I45.mat'));
I45=fI45.I45;
fI135=matfile(strcat(data_source,'I135.mat'));
I135=fI135.I135;


% M�thode des moindres carr�s sans refitting pour reconstruction des vecteurs de Stokes
% à partir de l'image polarimétrique de synthèse bruitée

M=size(I0,1);
N=size(I0,2);
MN = M*N;

% variables primales
S_hat=zeros(M,N,3);
S0_hat=zeros([1 MN]);
S1_hat=zeros([1 MN]);
S2_hat=zeros([1 MN]);


% Mesure du temps de calcul
tic

% Calcul de A_tilde, pseudo-inverse de A
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
A_tilde = (A' * A)^-1 * A';

% reshape de I
I = [I0(:)';I90(:)';I45(:)';I135(:)'];

sol = zeros(1,3);
for i=1:MN
    sol = A_tilde * I(:,i);
    S0_hat(i) = sol(1);
    S1_hat(i) = sol(2);
    S2_hat(i) = sol(3);
end
S_hat(:,:,1)=reshape(S0_hat,[M,N]);
S_hat(:,:,2)=reshape(S1_hat,[M,N]);
S_hat(:,:,3)=reshape(S2_hat,[M,N]);

% ########### REFITTING ############

% variables intermédiaires z_hat, utiles pour le co-support
% On a pos� p_hat = 0 et sigma_bi = gamma_bi = 1.
z_hat=gradient(S_hat);

% récupération du co-support
norm1=vecnorm(z_hat(:,:,1:2),2,3); % matrice des normes R2 sur le canal z_0 pixel-wise
nsup11=norm1>1; % indices où norme > 1 matrice 2D
nsup1r1=repmat(nsup11,1,1,2); % indices où norme > 1 matrice 3D
norm2=vecnorm(z_hat(:,:,3:6),2,3); % matrice des normes R4 sur les canaux (z_1,z_2) pixel-wise
nsup12=norm2>1; % indices où norme > 1 matrice 2D
nsup1r2=repmat(nsup12,1,1,6); % indices où norme > 1 matrice 3D
nsup1r2(:,:,1:2)=0;

% variables primales, problème de refitting
S_til=zeros(M,N,3);
S_new_til=zeros(M,N,3);
S_bar_til=S_til;

% variables duales, problème de refitting
p_til=zeros(M,N,6);

% paramètre b pour le prox de S_til
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
b=A'*[I0(:)';I90(:)';I45(:)';I135(:)'];

for k=1:Nb_iter_max
    
    % Exécution de l'itération, on utilise le mod�le de refitting avec TV
    % coupl�e sur les canaux (z_1,Z_2)
    [p_til,S_new_til]=cm_iteration_re(p_til,S_til,S_bar_til,sigma_re,nsup1r1,nsup1r2,z_hat,sigma_bi,tau_re,gamma_re,mu_re,M,N,MN,b);
    
    % Accélération de la convergence en continuant dans la direction
    % donnée par la mise à jour de S_til  
    S_bar_til=S_new_til+theta_re*(S_new_til-S_til);
   
    S_til=S_new_til;
    
end

% Mesure du temps de calcul
toc


% Sauvegarde des variables
save(strcat(data_dest,'S_hat.mat'),'S_hat');
save(strcat(data_dest,'S_til_mu_re_',sprintf('%.1f',mu_re),'.mat'),'S_til');
