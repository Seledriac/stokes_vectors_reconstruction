
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

% paramètres A et b pour le prox de S_hat
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
b=A'*[I0(:)';I90(:)';I45(:)';I135(:)'];


% Mesure du temps de calcul
tic

for k=1:Nb_iter_max
    
    % Ex�cution de l'it�ration
    S_new_hat=gma_iteration_bi(S_hat,M,N,MN,tau_bi,mu_bi,b);
    
    % Accélération de la convergence en continuant dans la direction
    % donnée par la mise à jour de S_hat
    S_bar_hat=S_new_hat+theta_bi*(S_new_hat-S_hat);    
   
    S_hat=S_new_hat;
    
end

% Mesure du temps de calcul
toc


% Sauvegarde des variables
save(strcat(data_dest,'S_hat_mu_bi_',sprintf('%.1f',mu_bi),'.mat'),'S_hat');
