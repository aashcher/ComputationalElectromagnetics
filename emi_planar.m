%{
Copyright © 2020 Alexey A. Shcherbakov. All rights reserved.

This file is part of ComputationalElectrodynamics.

ComputationalElectrodynamics is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

ComputationalElectrodynamics is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ComputationalElectrodynamics. If not, see <https://www.gnu.org/licenses/>.
%}
%% description:
% analyze the spetral and angular emission of a simple OLED structure
% on the basis of the S-matrix method
%% implementation
% STRUCTURE considered (simple green Alq3 OLED):
%   Air, semi-infinite
%   Al cathode, 150 nm
%   LiF, 0.8 nm
%   Alq3 green emitting layer, 50 nm
%   PEDOT:PSS, 80 nm
%   ITO anode, 90 nm
%   SiO2 substrate, semi-infinite
clc;
format long;

	% thicknesses of structure layers in micrometers:
l_thickness = [0.15, 0.0008, 0.05, 0.08, 0.09];
h_src = 0.025; % position of the emitting dipole layer
n_src = 3; % index of the emitting layer (Alq3)
% note: h_src sould be less than l_thickness(n_src)
	% read n-k data and normalize wavelengths where needed:
nk_Al = dlmread("./data/Al.nk");
nk_LiF = dlmread("./data/LiF.nk");
nk_Alq3 = dlmread("./data/Alq3.nk");
n_PP = dlmread("./data/PEDOT_PSS_n.nk");
k_PP = dlmread("./data/PEDOT_PSS_k.nk");
nk_ITO = dlmread("./data/ITO.nk");
nk_SiO2 = dlmread("./data/SiO2.nk");

nk_Al(:,1) = 0.001*nk_Al(:,1);
nk_LiF(:,1) = 0.001*nk_LiF(:,1);
nk_Alq3(:,1) = 0.001*nk_Alq3(:,1);
nk_ITO(:,1) = 0.001*nk_ITO(:,1);
nk_SiO2(:,1) = 0.001*nk_SiO2(:,1);

	% read and normalize Alq3 spectrum:
emi_Alq3 = dlmread("./data/Alq3.emi");
emi_Alq3(:,1) = 0.001*emi_Alq3(:,1);

	% parameters
wavelengths = 0.4:0.001:0.8;
k_vectors = 0.0001:0.01:1.5;

data_spec = [];

	% wavelength and k-space loop:
for i_wl = 1:numel(wavelengths)
	wl = wavelengths(i_wl);
	disp(wl);
		% normalized thicknesses:
	kh = l_thickness * (2*pi/wl);
	kh_src = h_src * (2*pi/wl);
		% calculate interpolated permittivities for the current wavelength:
	eps_Al = ( interp1(nk_Al(:,1), nk_Al(:,2), wl) + 1i * interp1(nk_Al(:,1), nk_Al(:,3), wl) )^2;
	eps_LiF = ( interp1(nk_LiF(:,1), nk_LiF(:,2), wl) + 1i * interp1(nk_LiF(:,1), nk_LiF(:,3), wl) )^2;
	eps_Alq3 = ( interp1(nk_Alq3(:,1), nk_Alq3(:,2), wl) + 1i * interp1(nk_Alq3(:,1), nk_Alq3(:,3), wl) )^2;
	eps_PP = ( interp1(n_PP(:,1), n_PP(:,2), wl) + 1i * interp1(k_PP(:,1), k_PP(:,2), wl) )^2;
	eps_ITO = ( interp1(nk_ITO(:,1), nk_ITO(:,2), wl) + 1i * interp1(nk_ITO(:,1), nk_ITO(:,3), wl) )^2;
	eps_SiO2 = ( interp1(nk_SiO2(:,1), nk_SiO2(:,2), wl) + 1i * interp1(nk_SiO2(:,1), nk_SiO2(:,3), wl) )^2;
	eps = [1, eps_Al, eps_LiF, eps_Alq3, eps_PP, eps_ITO, eps_SiO2];
		% source specral weight:
	weight = interp1(emi_Alq3(:,1), emi_Alq3(:,2), wl, 'linear', 0);
		% loop over k-space
	pow_sub = 0;
	pow_air = 0;
	for i_k = 1:numel(k_vectors)
		kx = k_vectors(i_k);
			% TE and TM scattering matrices below and above the emitting layer:
		[SM1e, SM1h, SM2e, SM2h] = get_src_smatrices(kx, kh, eps, kh_src, n_src);
			% power emitted to the substrate and to the air behind the substrate:
		if kx < real(sqrt(eps(end)))
			[P_sub, P_air] = get_power(kx, eps(n_src+1), eps(end), weight, SM1e, SM1h, SM2e, SM2h);
			pow_sub = pow_sub + P_sub;
			if kx < 1
				pow_air = pow_air + P_air;
			end
		end
	end
		% store data
	data_spec = [data_spec; [wl, pow_sub, pow_air]];
end

plot(data_spec(:,1), data_spec(:,2), "-b", ...
		 data_spec(:,1), data_spec(:,3), "-r", ...
		 emi_Alq3(:,1), 150*emi_Alq3(:,2), "-g");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function get_power
% calculate the power emitted by a dipole layer placed inside a layered
% structure
% input:
%		kx: in-plane wavevector projection normalized by the vacuum wavenumber
%		eps: permittivity of the emitting layer
%		sps_sub: substrate permittivity
%		weight: source spectral weigth
%		SM1e, SM1h: 2-by-2 S-matrices for the sub-structure below the dipole
%			layer for the TE and TM polarizations (NOTE: the structure is turned over)
%		SM2e, SM2h: 2-by-2 S-matrices for the sub-structure above the dipole
%			layer for the TE and TM polarizations (NOTE: the structure is turned over)
% output:
%		P_sub, P_air: power emitted to the substrate and to the air under
%			assumption that the reflections inside the substrate are non-coherent
function [P_sub, P_air] = get_power(kx, eps, eps_sub, weight, SM1e, SM1h, SM2e, SM2h)
	kz = sqrt(eps - kx^2);
	if angle(kz) < -1e-12; kz = -kz; end
		% dipole amplitudes:
	c_ep = 1/kz;
	c_em = 1/kz;
	c_h1p = 1;
	c_h1m = -1;
	c_h2p = kx/kz;
	c_h2m = kx/kz;
		% amplitudes of the field emitted to the substrate:
	ae = SM2e(2,1)*( c_ep + SM1e(2,2)*c_em ) / (1 - SM1e(2,2)*SM2e(1,1));
	ah1 = SM2h(2,1)*( c_h1p + SM1h(2,2)*c_h1m ) / (1 - SM1h(2,2)*SM2h(1,1));
	ah2 = SM2h(2,1)*( c_h2p + SM1h(2,2)*c_h2m ) / (1 - SM1h(2,2)*SM2h(1,1));
		% power carried by the plane waves emitted to the substrate:
	kz_sub = sqrt(eps_sub - kx^2);
	if angle(kz_sub) < -1e-12; kz_sub = -kz_sub; end
	P_sub = 0.5*weight*( norm(ae)*real(kz_sub) + (norm(ah1) + norm(ah2))*real(kz_sub/eps_sub) );
		% S-matrices of the whole structure:
	SMe = smatrix_multiply(SM1e, SM2e);
	SMh = smatrix_multiply(SM1h, SM2h);
		% S-matrices of the substrate-air interface:
	Se = smatrix_interface(eps_sub, 1, kx, "TE");
	Sh = smatrix_interface(eps_sub, 1, kx, "TM");
		% power reflection coefficients at the structure intefrace and at the
		% substrate-air interface:
	RSe = norm(SMe(2,2));
	RSh = norm(SMh(2,2));
	RAe = norm(Se(1,1));
	RAh = norm(Sh(1,1));
		% power transmssion coefficients for the emitted radiation:
	if kx < 1
		N = 4*real(kz_sub)/kx; % estimate of the number of reflections inside the substrate
		Ce = (1 - RAe) * (1 - exp(N*log(RAe*RSe))) / (1 - RAe*RSe);
		Ch = (1 - RAh) * (1 - exp(N*log(RAh*RSh))) / (1 - RAh*RSh);
		P_air = 0.5*weight*( Ce*norm(ae)*real(kz_sub) + Ch*(norm(ah1) + norm(ah2))*real(kz_sub/eps_sub) );
	else
		P_air = 0;
	end
end

%% function get_src_smatrices
% function for calculation of upper and lower S-matrices relative to a
% dipole source located inside a layered structure (NOTE: the structure is turned over)
% input:
%		kx: in-plane wavevector projection normalized by the vacuum wavenumber
%		kh: array of normalized thicknesses of structure layers
%		eps: array of permittivities for the superstrate, structure layers, and
%			the substrate
%		kh_src: depth of the dipole layer inside the emitting layer
%		n_src: index of the emitting layer within the structure (corresponds to
%			kh(n_src) and eps(n_src+1))
% output:
%		SM1e, SM1h: 2-by-2 S-matrices for the sub-structure below the dipole
%			layer for the TE and TM polarizations
%		SM2e, SM2h: 2-by-2 S-matrices for the sub-structure above the dipole
%			layer for the TE and TM polarizations
function [SM1e, SM1h, SM2e, SM2h] = get_src_smatrices(kx, kh, eps, kh_src, n_src)
	if n_src > numel(kh)
		error("get_smatrices: number of structure layers is less than the source layer number");
	end
	if numel(eps) - numel(kh) ~= 2
		error("get_smatrices: incompatible size of input arrays");
	end
		% S-matrices below the dipole layer:
	SM1e = smatrix_interface(eps(1), eps(2), kx, "TE");
	SM1h = smatrix_interface(eps(1), eps(2), kx, "TM");
	for i = 1:n_src-1
		SM_tmp = smatrix_layer(kh(i), eps(i+1), kx);
		SM1e = smatrix_multiply(SM1e, SM_tmp);
		SM1h = smatrix_multiply(SM1h, SM_tmp);
		SM_tmp = smatrix_interface(eps(i+1), eps(i+2), kx, "TE");
		SM1e = smatrix_multiply(SM1e, SM_tmp);
		SM_tmp = smatrix_interface(eps(i+1), eps(i+2), kx, "TM");
		SM1h = smatrix_multiply(SM1h, SM_tmp);
	end
	SM_tmp = smatrix_layer(kh_src, eps(n_src+1), kx);
	SM1e = smatrix_multiply(SM1e, SM_tmp);
	SM1h = smatrix_multiply(SM1h, SM_tmp);
		% S-matrices above the dipole layer:
	SM2e = smatrix_layer(kh(n_src)-kh_src, eps(n_src+1), kx);
	SM2h = SM2e;
	SM_tmp = smatrix_interface(eps(n_src+1), eps(n_src+2), kx, "TE");
	SM2e = smatrix_multiply(SM2e, SM_tmp);
	SM_tmp = smatrix_interface(eps(n_src+1), eps(n_src+2), kx, "TM");
	SM2h = smatrix_multiply(SM2h, SM_tmp);
	for i = n_src+1:numel(kh)
		SM_tmp = smatrix_layer(kh(i), eps(i+1), kx);
		SM2e = smatrix_multiply(SM2e, SM_tmp);
		SM2h = smatrix_multiply(SM2h, SM_tmp);
		SM_tmp = smatrix_interface(eps(i+1), eps(i+2), kx, "TE");
		SM2e = smatrix_multiply(SM2e, SM_tmp);
		SM_tmp = smatrix_interface(eps(i+1), eps(i+2), kx, "TM");
		SM2h = smatrix_multiply(SM2h, SM_tmp);
	end
end

%%% END %%%