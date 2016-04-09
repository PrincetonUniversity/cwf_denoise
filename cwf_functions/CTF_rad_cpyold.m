
function [ctf_rad] = CTF_rad_cpyold(n, res, lambda, defocus, Cs, B, alpha, r)
% f = ctf(n, res, lambda, defocus, Cs, B, alpha)
% Compute the contrast transfer function corresponding to the resolution
% res in A/pixel.  Lambda is in A, Cs is in mm,
% B in A^2 and alpha in radians.  Defocus is in microns.
% The result is returned in an nxn matrix with h(n/2+1) corresponding
% to the zero-frequency amplitude.  Thus you must use fftshift() on
% the result before performing a Fourier transform.  For example,
% to simulate the effect of the CTF on an image m, do this:
% fm=cfft2(m);
% cm=icfft2(fm.*h));

% Cs term fixed.  fs 4 Apr 04
%
% The first zero occurs at lambda*defocus*f0^2=1.
% e.g. when lambda=.025A, defocus=1um, then f0=1/16�.

% Usage: lambda is the wavelength determined by the voltage. For example,
% lambda = EWavelength (300) computes the wavelength in A for 300KV.
% n is the size of the image, for example n=129. 
% res = 3 A/pixel.
% defocus is in micrometers, usually it is in the range of 1 to 4 micrometers.
% Cs is normally 2.
% B determines the decay envelope and it can be 0, 10 or 100.
% alpha is usually set at 0.07.


f0 = 1/(n*res);  % Spatial frequency unit (inverse �)

    k2=-pi*lambda*defocus*1e4*f0^2;
    k4= pi/2*Cs*lambda^3*1e7*f0^4;  %2.19^4 in Frank's book

kr=f0^2*B;  % B-factor

r2=r.^2;
if Cs==0
    ctf_rad=sin(k2*r2-alpha).*exp(-kr*r2);
else
    ctf_rad=sin(k2*r2+k4*r2.*r2-alpha).*exp(-kr*r2);
end
