
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
    sol = A_tilde * I(i);
    S0_hat(i) = sol(1);
    S1_hat(i) = sol(2);
    S2_hat(i) = sol(3);
end
S_hat(:,:,1)=reshape(S0_hat,[M,N]);
S_hat(:,:,2)=reshape(S1_hat,[M,N]);
S_hat(:,:,3)=reshape(S2_hat,[M,N]);

% Mesure du temps de calcul
toc


% Sauvegarde des variables
save(strcat(data_dest,'S_hat.mat'),'S_hat');


