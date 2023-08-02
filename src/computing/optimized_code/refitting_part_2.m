% Renvoie le point proximal de z_til par rapport à la fonction 
% phi_SD^*(.,psi,I_hat) pour un pixel appartenant au co-support.
% Le résultat est la projection de z_til sur la boule de centre
% -psi/norm(psi) de rayon gamme_re (en norme euclidienne)
% (la dimension de la boule dépend des dimensions des paramètres).
% dimension fixée par le paramètre dim

function sol_ex = refitting_part_2(z_til,psi,gamma_re,dim)
norm_psi = repmat(vecnorm(psi,2,1),dim,1,1);
psi_normalized = psi ./ norm_psi;
tmp_sum = z_til + gamma_re * psi_normalized;
norm_tmp_sum = repmat(vecnorm(tmp_sum,2,1),dim,1,1);
tmp_sum(norm_tmp_sum>gamma_re) = gamma_re * tmp_sum(norm_tmp_sum>gamma_re) ./ norm_tmp_sum(norm_tmp_sum>gamma_re);
sol_ex = tmp_sum - gamma_re * psi_normalized;