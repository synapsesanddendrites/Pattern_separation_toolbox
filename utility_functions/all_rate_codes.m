% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [all_codes] = all_rate_codes(max_value, n_trial)

basis_vec=0:max_value;
if n_trial==1
    all_codes=combvec(basis_vec)';
elseif n_trial==2
    all_codes=combvec(basis_vec,basis_vec)';
elseif n_trial==3
    all_codes=combvec(basis_vec,basis_vec,basis_vec)';
elseif n_trial==4
    all_codes=combvec(basis_vec,basis_vec,basis_vec,basis_vec)';
elseif n_trial==5
    all_codes=combvec(basis_vec,basis_vec,basis_vec,basis_vec,basis_vec)';
elseif n_trial==6
    all_codes=combvec(basis_vec,basis_vec,basis_vec,basis_vec,basis_vec,basis_vec)';
elseif n_trial==7
    all_codes=combvec(basis_vec,basis_vec,basis_vec,basis_vec,basis_vec,basis_vec,basis_vec)';
elseif n_trial==8
    all_codes=combvec(basis_vec,basis_vec,basis_vec,basis_vec,basis_vec,basis_vec,basis_vec,basis_vec)';
end

end