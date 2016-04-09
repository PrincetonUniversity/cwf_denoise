n=256;
K=1000;
addpath ~/Dropbox/asinger/80S/80s_Tejal/SUGAR/sugar-master/my_scripts/ccwf_v6/final_scripts/expt_ip3
[projections] = emd5278_proj_full(n,K);
save /scratch/tbhamre/expt_ip3/ip3_proj.mat  projections
