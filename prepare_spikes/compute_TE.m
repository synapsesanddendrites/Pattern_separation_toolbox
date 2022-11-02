% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [TE_sum] = compute_TE(joint_distribution,delay_distribution) 
dd_norm=sum(delay_distribution(:));
jd_norm=sum(joint_distribution(:));

sz=size(joint_distribution);
H1=0;
for lag_ind=1:sz(3)
    PZ=delay_distribution(lag_ind)/dd_norm;
    for post_ind=1:sz(2)
        PYZ=sum(joint_distribution(:,post_ind,lag_ind),'all')/jd_norm;
        if PYZ>0
            H1=H1-PYZ*log2(PYZ/PZ);
        end
    end
end

H2=0;
for lag_ind=1:sz(3)
    for pre_ind=1:sz(1)
        PXZ=sum(joint_distribution(pre_ind,:,lag_ind),'all')/jd_norm;
        for post_ind=1:sz(2)
            PXYZ=joint_distribution(pre_ind,post_ind,lag_ind)/jd_norm;
            if PXYZ>0
                H2=H2-PXYZ*log2(PXYZ/PXZ);
            end
        end
    end
end
TE_sum=H1-H2;
end