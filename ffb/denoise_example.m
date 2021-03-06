%Denoise example code, input image
var_hat = 1; %noise variance, change accordingly
n_im = 3; % number of denoised images
[ c, R ] = avg_pspec(data, var_hat); %Estimate band limit and compact support size
[ timing, coeff, mean_coeff, sPCA_coeff, U, D ] = jobscript_FFBsPCA(data, c, R, var_hat); % Compute Fast steerable PCA.
[ fn ] = IFT_FB(R, c); % compute inverse fourier transform of the Fourier-Bessel functions
L0=size(data,1);
[ mean_image, denoised ] = denoise_images_analytical(U, fn, mean_coeff, sPCA_coeff, L0, R, n_im); %generate denoised images.


