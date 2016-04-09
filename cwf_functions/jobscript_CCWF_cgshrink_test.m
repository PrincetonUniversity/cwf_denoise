function [ denoised_coeff ] = jobscript_CCWF_cgshrink_test(index, CTF_rad_all, basis, sample_points,  mean_coeff, coeff_pos_k,  noise_v)

% Main function to estimate covariance (white noise) and obtain denoised image coefficients by deconvolution (Wiener filter)
% Inputs:
% index: CTF indices for images
% w_f_rad: Whitening filter evaluated along radial quadrature points
% CTF_rad_all: matrix of all radial CTF's, each column is a distinct CTF
% basis: precomputed Fourier Bessel basis
% sample_points: Quadrature points
% mean_coeff: FB coefficients of the mean image
% coeff_pos_k: FB coefficients of images for non-negative angular frequencies
% noise_v: Noise variance
% OUTPUTS:
% denoised_coeff: Cell of FB coeffs of denoised images, each cell corresponds to a non-negative angular frequency
% eigims: Cell of eigenimages corresponding to each angular frequency
% egval: EIgenvalues corresponding to each angular frequency block
% num_eigs: Number of non-zero eigenvalues for each block
% Tejal Jan 2016

nim=length(index);
ang_freqs=basis.ang_freqs;

for i=1:max(index)
    weight(i)=length(find(index==i));
end

denoised_coeff=cell(length(coeff_pos_k), 1);
eigims=cell(length(coeff_pos_k), 1);
% Only pos k, to be compatibel with Jane's ffb
for k=unique(ang_freqs)'
    
    if k>=0
        tmp=coeff_pos_k{k+1};
        %tr_c=tr_coeff_pos_k{k+1};
        %tr_cov_k=real(tr_c*tr_c')/(1e5);
        % tr_sample_k=tr_samplecov{k+1};
        %tmp1=conj(tmp);
        clear DD A
        for i=1:max(index)
            
            DD{i}=real(tmp(:,find(index==i))*tmp(:,find(index==i))'/length(find(index==i)));
            [UDD,SDD]=eig((DD{i}));
            
            % Shrinkage for sample covariance
            n_i=length(find(index==i));
            if (k==0)
                gamma=size(DD{i},1)/n_i;
            else
                gamma=size(DD{i},1)/(2*n_i);
            end
            
            
            % tr_DD{i}=real(tr_sample_k(:,find(index==i))*tr_sample_k(:,find(index==i))'/length(find(index==i)));
            
            A{i}=calc_fb_CTF(CTF_rad_all(:,i),basis.Phi_ns{k+1}, sample_points);
        end
        
        
        %No regularization for white noise case
        regu=0;
        
        [C]=solve_cg_shrink(A,DD,weight,noise_v, nim,k,regu);  % Changed, solve_LS_c from solve_LS which doesn't include weights
        
        
        %% Ensure positivity of the estimated covariance matrix
        [U,S]=eig(C);
        C=U*diag(max(diag(S),0))*U';
        
        for i=1:max(index)
            
            h_posk = C*A{i}'*inv(A{i}*C*A{i}'+noise_v*eye(size(A{i},1)));
            
            if (k==0)
                denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) + repmat(mean_coeff{k+1},1,numel(find(index==i))); %disp('Subtracting mean in deconvolution')
            else
                denoised_coeff{k+1}(:,find(index==i)) = h_posk*tmp(:,find(index==i)) ; %disp('Subtracting mean in deconvolution')
            end
            
        end
    end
end
