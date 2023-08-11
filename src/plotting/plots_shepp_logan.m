close all;
% clc;


% flags
model_1 = 'cm';
mu_bi_1 = 7.0;
% mu_re_1 = 1.5;
fmodel_2 = false;
if(fmodel_2)
    model_2 = 'cm';
    mu_bi_2 = 1.0;
    mu_re_2 = 1.0;
end


% RÃ©cupÃ©ration de l'image polarimétrique
% I clean
fI1=matfile('../../data_in/shepp_logan_phantom/I0_clean.mat');
I1=fI1.I1;
fI2=matfile('../../data_in/shepp_logan_phantom/I90_clean.mat');
I2=fI2.I2;
fI3=matfile('../../data_in/shepp_logan_phantom/I45_clean.mat');
I3=fI3.I3;
fI4=matfile('../../data_in/shepp_logan_phantom/I135_clean.mat');
I4=fI4.I4;
% I bruité
fI0=matfile('../../data_in/shepp_logan_phantom/I0.mat');
I0=fI0.I0;
fI90=matfile('../../data_in/shepp_logan_phantom/I90.mat');
I90=fI90.I90;
fI45=matfile('../../data_in/shepp_logan_phantom/I45.mat');
I45=fI45.I45;
fI135=matfile('../../data_in/shepp_logan_phantom/I135.mat');
I135=fI135.I135;
M=size(I0,1);
N=size(I0,2);
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];


% Récupération des vecteurs de Stokes
% Pour le modèle 1
if ( strcmp(model_1,'cm') )
    fS_hat=matfile(strcat('../../data_out/cm/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
    S_model_1=fS_hat.S_hat;
elseif ( strcmp(model_1,'cm_r') )
    fS_til=matfile(strcat('../../data_out/cm/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
    S_model_1=fS_til.S_til;
elseif ( strcmp(model_1,'dm') )
    fS_hat=matfile(strcat('../../data_out/dm/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
    S_model_1=fS_hat.S_hat;
elseif ( strcmp(model_1,'dm_r') )
    fS_til=matfile(strcat('../../data_out/dm/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
    S_model_1=fS_til.S_til;
elseif ( strcmp(model_1,'gma') )
    fS_hat=matfile(strcat('../../data_out/gma/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
    S_model_1=fS_hat.S_hat;
elseif ( strcmp(model_1,'gma_r') )
    fS_til=matfile(strcat('../../data_out/gma/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
    S_model_1=fS_til.S_til;
elseif ( strcmp(model_1,'ls') )
    fS_hat=matfile(strcat('../../data_out/ls/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
    S_model_1=fS_hat.S_hat;
elseif ( strcmp(model_1,'ls_r') )
    fS_til=matfile(strcat('../../data_out/ls/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
    S_model_1=fS_til.S_til;
end
S0_model_1=S_model_1(:,:,1);
S1_model_1=S_model_1(:,:,2);
S2_model_1=S_model_1(:,:,3);
% Pour le modèle 2
if(fmodel_2)
    if ( strcmp(model_2,'cm') )
        fS_hat=matfile(strcat('../../data_out/cm/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
        S_model_2=fS_hat.S_hat;
    elseif ( strcmp(model_2,'cm_r') )
        fS_til=matfile(strcat('../../data_out/cm/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
        S_model_2=fS_til.S_til;
    elseif ( strcmp(model_2,'dm') )
        fS_hat=matfile(strcat('../../data_out/dm/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
        S_model_2=fS_hat.S_hat;
    elseif ( strcmp(model_2,'dm_r') )
        fS_til=matfile(strcat('../../data_out/dm/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
        S_model_2=fS_til.S_til;
    elseif ( strcmp(model_2,'gma') )
        fS_hat=matfile(strcat('../../data_out/gma/shepp_logan_phantom/S_hat_mu_bi_',sprintf('%.1f',mu_bi_1),'.mat'));
        S_model_2=fS_hat.S_hat;
    elseif ( strcmp(model_2,'gma_r') )
        fS_til=matfile(strcat('../../data_out/gma/shepp_logan_phantom/S_til_mu_bi_',sprintf('%.1f',mu_bi_1),'_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
        S_model_2=fS_til.S_til;
    elseif ( strcmp(model_2,'ls') )
        fS_hat=matfile(strcat('../../data_out/ls/shepp_logan_phantom/S_hat.mat'));
        S_model_2=fS_hat.S_hat;
    elseif ( strcmp(model_2,'ls_r') )
        fS_til=matfile(strcat('../../data_out/ls/shepp_logan_phantom/S_til_mu_re_',sprintf('%.1f',mu_re_1),'.mat'));
        S_model_2=fS_til.S_til;
    end
    S0_model_2=S_model_2(:,:,1);
    S1_model_2=S_model_2(:,:,2);
    S2_model_2=S_model_2(:,:,3);
end


% Reconstruction de l'image polarimÃ©trique d'origine à partir des vecteurs
% de Stokes
% Modèle 1
I0_recovered_model_1 = zeros(M,N);
I90_recovered_model_1 = zeros(M,N);
I45_recovered_model_1 = zeros(M,N);
I135_recovered_model_1 = zeros(M,N);
for i=1:M
    for j=1:N
        I_recovered_model_1=A*[S0_model_1(i,j);S1_model_1(i,j);S2_model_1(i,j)];
        I0_recovered_model_1(i,j)=I_recovered_model_1(1);
        I90_recovered_model_1(i,j)=I_recovered_model_1(2);
        I45_recovered_model_1(i,j)=I_recovered_model_1(3);
        I135_recovered_model_1(i,j)=I_recovered_model_1(4);
    end
end
% Modèle 2
if(fmodel_2)
    I0_recovered_model_2 = zeros(M,N);
    I90_recovered_model_2 = zeros(M,N);
    I45_recovered_model_2 = zeros(M,N);
    I135_recovered_model_2 = zeros(M,N);
    for i=1:M
        for j=1:N
            I_recovered_model_2=A*[S0_model_2(i,j);S1_model_2(i,j);S2_model_2(i,j)];
            I0_recovered_model_2(i,j)=I_recovered_model_2(1);
            I90_recovered_model_2(i,j)=I_recovered_model_2(2);
            I45_recovered_model_2(i,j)=I_recovered_model_2(3);
            I135_recovered_model_2(i,j)=I_recovered_model_2(4);
        end
    end
end


if(fmodel_2)


    % figure 1, I0 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1 VS reconstruit par le modèle 2, 1x4
    c_min=min([min(min(I1)), min(min(I0)), min(min(I0_recovered_model_1)), min(min(I0_recovered_model_2))]);
    c_max=max([max(max(I1)), max(max(I0)), max(max(I0_recovered_model_1)), max(max(I0_recovered_model_2))]);
    t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I0_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I0_recovered_model_2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    set(gcf, 'Position', get(0, 'Screensize'));
%     sgtitle('I0 : original VS altered VS model 1 VS model 2');
    filename = fullfile('../../figures/',strcat(model_1,'_',model_2,'_shepp_logan_I0_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 2, I90 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1 VS reconstruit par le modèle 2, 1x4
    c_min=min([min(min(I2)), min(min(I90)), min(min(I90_recovered_model_1)), min(min(I90_recovered_model_2))]);
    c_max=max([max(max(I2)), max(max(I90)), max(max(I90_recovered_model_1)), max(max(I90_recovered_model_2))]);
    t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I90,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I90_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I90_recovered_model_2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    set(gcf, 'Position', get(0, 'Screensize'));
%     sgtitle('I90 : original VS altered VS model 1 VS model 2');
    filename = fullfile('../../figures/',strcat(model_1,'_',model_2,'_shepp_logan_I90_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 3, I45 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1 VS reconstruit par le modèle 2, 1x4
    c_min=min([min(min(I3)), min(min(I45)), min(min(I45_recovered_model_1)), min(min(I45_recovered_model_2))]);
    c_max=max([max(max(I3)), max(max(I45)), max(max(I45_recovered_model_1)), max(max(I45_recovered_model_2))]);
    t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I3,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I45,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I45_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I45_recovered_model_2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    set(gcf, 'Position', get(0, 'Screensize'));
%     sgtitle('I45 : original VS altered VS model 1 VS model 2');
    filename = fullfile('../../figures/',strcat(model_1,'_',model_2,'_shepp_logan_I45_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 4, I135 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1 VS reconstruit par le modèle 2, 1x4
    c_min=min([min(min(I4)), min(min(I135)), min(min(I135_recovered_model_1)), min(min(I135_recovered_model_2))]);
    c_max=max([max(max(I4)), max(max(I135)), max(max(I135_recovered_model_1)), max(max(I135_recovered_model_2))]);
    t = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I4,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I135,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I135_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I135_recovered_model_2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    set(gcf, 'Position', get(0, 'Screensize'));
%     sgtitle('I135 : original VS altered VS model 1 VS model 2');
    filename = fullfile('../../figures/',strcat(model_1,'_',model_2,'_shepp_logan_I135_comparison.png'));
    exportgraphics(t,filename);
    
    
    % % Changement de colormap
    % m=1000;
    % cm_inferno=inferno(m);
    % colormap(inferno);
    
    
    % figure 5, DOLP : modèle 1 VS modèle 2, 1x2
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    dolp_model_1 = sqrt(S1_model_1.^2+S2_model_1.^2)./S0_model_1;
    dolp_model_2 = sqrt(S1_model_2.^2+S2_model_2.^2)./S0_model_2;
    c_min = min([min(min(dolp_model_1)),min(min(dolp_model_2))]);
    c_max = max([max(max(dolp_model_1)),max(max(dolp_model_2))]);
    nexttile,imagesc(dolp_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imagesc(dolp_model_2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    set(gcf, 'Position', get(0, 'Screensize'));
%     sgtitle('DOLP : original VS model 1 VS model 2');
    filename = fullfile('../../figures/',strcat(model_1,'_',model_2,'_shepp_logan_dolp_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 6, AOLP : modèle 1 VS modèle 2, 1x2
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    aolp_model_1 = 0.5*atan2(S2_model_1,S1_model_1);
    aolp_model_2 = 0.5*atan2(S2_model_2,S1_model_2);
    c_min = min([min(min(aolp_model_1)),min(min(aolp_model_2))]);
    c_max = max([max(max(aolp_model_1)),max(max(aolp_model_2))]);
    nexttile,imagesc(aolp_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imagesc(aolp_model_2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    set(gcf, 'Position', get(0, 'Screensize'));
%     sgtitle('AOLP : original VS model 1 VS model 2');
    filename = fullfile('../../figures/',strcat(model_1,'_',model_2,'_shepp_logan_aolp_comparison.png'));
    exportgraphics(t,filename);


else


    % figure 1, I0 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1, 1x3
    c_min=min([min(min(I1)), min(min(I0)), min(min(I0_recovered_model_1))]);
    c_max=max([max(max(I1)), max(max(I0)), max(max(I0_recovered_model_1))]);
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I0_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
%     sgtitle('I0 : original VS altered VS model 1');
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model_1,'_shepp_logan_I0_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 2, I90 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1, 1x3
    c_min=min([min(min(I2)), min(min(I90)), min(min(I90_recovered_model_1))]);
    c_max=max([max(max(I2)), max(max(I90)), max(max(I90_recovered_model_1))]);
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I90,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I90_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
%     sgtitle('I90 : original VS altered VS model 1');
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model_1,'_shepp_logan_I90_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 3, I45 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1, 1x3
    c_min=min([min(min(I3)), min(min(I45)), min(min(I45_recovered_model_1))]);
    c_max=max([max(max(I3)), max(max(I45)), max(max(I45_recovered_model_1))]);
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I3,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I45,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I45_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
%     sgtitle('I45 : original VS altered VS model 1');
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model_1,'_shepp_logan_I45_comparison.png'));
    exportgraphics(t,filename);
    
    
    % figure 4, I135 : d'origine VS bruitÃ©e VS reconstruit
    % par le modèle 1, 1x3
    c_min=min([min(min(I4)), min(min(I135)), min(min(I135_recovered_model_1))]);
    c_max=max([max(max(I4)), max(max(I135)), max(max(I135_recovered_model_1))]);
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    nexttile,imshow(I4,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I135,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
    nexttile,imshow(I135_recovered_model_1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
%     sgtitle('I135 : original VS altered VS model 1');
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model_1,'_shepp_logan_I135_comparison.png'));
    exportgraphics(t,filename);
    
    
    % % Changement de colormap
    % m=1000;
    % cm_inferno=inferno(m);
    % colormap(inferno);
    
    
    % figure 5, DOLP : modèle 1, 1x1
    f=figure();
    dolp_model_1 = sqrt(S1_model_1.^2+S2_model_1.^2)./S0_model_1;
    imagesc(dolp_model_1);colorbar;axis off;
%     sgtitle('DOLP : model 1');
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model_1,'_shepp_logan_dolp.png'));
    exportgraphics(f,filename);
    
    
    % figure 6, AOLP : modèle 1, 1x1
    f=figure();
    aolp_model_1 = 0.5*atan2(S2_model_1,S1_model_1);
    imagesc(aolp_model_1);colorbar;axis off;
%     sgtitle('AOLP : model 1');
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model_1,'_shepp_logan_aolp.png'));
    exportgraphics(f,filename);


end
