% Renvoie le point proximal de z_til par rapport à la fonction 
% phi_SD^*(.,psi,I_hat) pour un pixel appartenant au co-support.
% Le résultat est la projection de z_til sur la boule de centre
% -psi/norm(psi) de rayon gamme_re (en norme euclidienne)
% (la dimension de la boule dépend des dimensions des paramètres).

function sol_ex = refitting_part(z_til,psi,gamma_re)

norm_psi = norm(psi,2);
tmp_sum = z_til + gamma_re * (psi / norm_psi);
norm_tmp_sum = norm(tmp_sum,2);
sol_ex = gamma_re * ( ( tmp_sum / max( gamma_re, norm_tmp_sum ) ) - ( psi / norm_psi ) );