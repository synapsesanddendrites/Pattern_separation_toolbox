% the pattern separation toolbox: analyse spike train ensembles
% Copyright (C) 2022  Alexander D Bird

function [RR_sum] = compute_RR(joint_distribution,input_distribution,output_distribution)    
id_norm=sum(input_distribution(:));
od_norm=sum(output_distribution(:));
jd_norm=sum(joint_distribution(:));

PXY=joint_distribution/jd_norm;
PX=input_distribution/id_norm;
PY=output_distribution/od_norm;

[pre,post]=find(PXY);
linearidx = sub2ind(size(PXY),pre,post);

RR_sum=sum(PXY(linearidx).*log2(PXY(linearidx)./(PX(pre).*PY(post))),'all');
end
