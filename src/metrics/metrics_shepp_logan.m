
% Solution biaisée du mod�le appliqué à
% l'image polarimétrique d'entrée I
fS_hat=matfile(strcat('../../data_out/',model,'/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi),'.mat'));
S_hat=fS_hat.S_hat;
S0=S_hat(:,:,1);
S1=S_hat(:,:,2);
S2=S_hat(:,:,3);

if refitting == true
    % Solution débiaisée du mod�le appliqué à
    % l'image polarimétrique d'entrée I
    fS_til=matfile(strcat('../../data_out/',model,'/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi),'_mu_re_',sprintf('%.1f',mu_re),'.mat'));
    S_til=fS_til.S_til;
    S0_refitted=S_til(:,:,1);
    S1_refitted=S_til(:,:,2);
    S2_refitted=S_til(:,:,3);
end
    
% Reconstruction de l'image polarimétrique à partir du
% S reconstruit biaisé
for i=1:M
    for j=1:N
        I_recovered=A*[S0(i,j);S1(i,j);S2(i,j)];
        I0_recovered(i,j)=I_recovered(1);
        I90_recovered(i,j)=I_recovered(2);
        I45_recovered(i,j)=I_recovered(3);
        I135_recovered(i,j)=I_recovered(4);
    end
end
% PSNR entre image polarimétrique non bruitée et reconstruite biaisée
% L'échelle est de 0 à 1 pour chaque canal et il y a 4MN pixels (car MN
% superpixels)
PSNR=10*log10(1/((1/(4*M*N))*(sum(sum((I1-I0_recovered).^2))+sum(sum((I2-I90_recovered).^2))+...
    sum(sum((I3-I45_recovered).^2))+sum(sum((I4-I135_recovered).^2)))));
% SSIM entre image polarimétrique non bruitée et reconstruite biaisée
SSIM=[ssim(I1,I0_recovered);
ssim(I2,I90_recovered);
ssim(I3,I45_recovered);
ssim(I4,I135_recovered)
];

if refitting == true
    % Reconstruction de l'image polarimétrique à partir du
    % S reconstruit débiaisé
    for i=1:M
        for j=1:N
            I_recovered_refitted=A*[S0_refitted(i,j);S1_refitted(i,j);S2_refitted(i,j)];
            I0_recovered_refitted(i,j)=I_recovered_refitted(1);
            I90_recovered_refitted(i,j)=I_recovered_refitted(2);
            I45_recovered_refitted(i,j)=I_recovered_refitted(3);
            I135_recovered_refitted(i,j)=I_recovered_refitted(4);
        end
    end
    % PSNR entre image polarimétrique non bruitée et reconstruite débiaisée
    % L'échelle est de 0 à 1 pour chaque canal et il y a 4MN pixels (car MN
    % superpixels)
    PSNR_refitted=10*log10(1/((1/(4*M*N))*(sum(sum((I1-I0_recovered_refitted).^2))+sum(sum((I2-I90_recovered_refitted).^2))+...
        sum(sum((I3-I45_recovered_refitted).^2))+sum(sum((I4-I135_recovered_refitted).^2)))));
    % SSIM entre image polarimétrique non bruitée et reconstruite débiaisée
    SSIM_refitted=[ssim(I1,I0_recovered_refitted);
    ssim(I2,I90_recovered_refitted);
    ssim(I3,I45_recovered_refitted);
    ssim(I4,I135_recovered_refitted)
    ];
end

% note : SSIM incorpore dans sa mesure :
% - la différence de luminance, qui correspond à la différence de moyenne d'intensités entre Ix et Ix_recovered
% - la différence de contraste, qui correspond à la différence de variance
% - la différence de similarité de structure, qui correspond au coefficient
% de corrélation entre I_x et I_recovered
