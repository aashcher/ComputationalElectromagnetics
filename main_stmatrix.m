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
% calculate reflection spectra from inhomogeneous slabs; from periodically
% layered slabs; from homogeneous anisotropic slabs
%% implementation:
clc;
clear;
%% initialization
wl = 0.55; % wavelength
wv = 2*pi/wl; % vacuum wavenumber
theta = pi/9; % propagation angle in the air
kx0 = sin(theta); % in-plane wavevector projection
pol = 'TE'; % polarization, either 'TE' or 'TM'

%% reflection from a slab with sinusoidally varying refractive index
% {
n_sl = 10; % number of slices
H = 0.25; % slab thickness
dh = H/n_sl;

wls = 0.4:0.001:0.8; % row of wavelengths for spectrum simulation
reflection_ST = zeros(numel(wls),2);

for iwl = 1:numel(wls)
	kdh = 2*pi*dh / wls(iwl);
	SM = zeros(2,2); SM(1,2) = 1; SM(2,1) = 1; % S-matrix of nothing
	TM = zeros(2,2); TM(1,1) = 1; TM(2,2) = 1; % T-matrix of nothing
	eps_2 = 1;
		% accumulate S-matrix of the slab
	for is = 1:n_sl
		ri = 1 + 10*sin(pi*(is+0.5)/n_sl); % sinusoidally varying refractive index
		eps_1 = eps_2;
		eps_2 = ri^2;
			% calculate and accumulate slab S- and T- matrices
		SM = smatrix_multiply(SM, smatrix_interface(eps_1,eps_2,kx0,pol));
		SM = smatrix_multiply(SM, smatrix_layer(kdh,eps_2,kx0));
		TM = TM * tmatrix_interface(eps_1,eps_2,kx0,pol);
		TM = TM * tmatrix_layer(kdh,eps_2,kx0);
	end
		% the last interface
	ri = 1;
	eps_1 = eps_2;
	eps_2 = ri^2;
	SM = smatrix_multiply(SM, smatrix_interface(eps_1,eps_2,kx0,pol));
	TM = TM * tmatrix_interface(eps_1,eps_2,kx0,pol);
		% store data
	reflection_ST(iwl,1) = abs(SM(1,1)^2);
	reflection_ST(iwl,2) = abs((TM(2,1)/TM(2,2))^2); % reflection coefficient from T-matrix
end
	% plot reflection spectrum
plot(wls,reflection_ST(:,1));
hold on;
plot(wls,reflection_ST(:,2));
hold off;
xlabel('\lambda');
ylabel('R');

return;
%}
%% reflection from a periodic multilayer
%{
kx0 = 0.6; % in-plane wavevector projection
pol = 'TE'; % polarization, either 'TE' or 'TM'

t1 = 0.2; % thickness of the first layer in length units
t2 = 0.1; % thickness of the second layer in length units
ri1 = 4 + 0*1i; % refractive index of the first layer
ri2 = 1 + 0*1i; % refractive index of the second layer
	% permittivities:
eps1 = ri1^2;
eps2 = ri2^2;
eps = 1; % surrounding medium

np = 10; % number of periods
H = np*(t1 + t2); % total slab thickness

  % effective medium constants
epsR_o = (t1*eps1 + t2*eps2)/(t1 + t2);
epsR_e = (t1 + t2)/(t1/eps1 + t2/eps2);

wls = 0.4:0.0005:0.8; % row of wavelengths for spectrum simulation in length units
reflection_ST = zeros(numel(wls),3);

	% spectrum loop
for iwl = 1:numel(wls)
	wv = 2*pi/wls(iwl); % vacuum wavenumber
	%eps2 = get_epsAu_Drude(wls(iwl));
		% normalized thicknesses:
	kh1 = wv*t1;
	kh2 = wv*t2;
	kh = (kh1 + kh2)*np;

		% S- and T- matrices of one single period
	SMP = smatrix_interface(eps2,eps1,kx0,pol);
	SMP = smatrix_multiply(SMP, smatrix_layer(kh1,eps1,kx0));
	SMP = smatrix_multiply(SMP, smatrix_interface(eps1,eps2,kx0,pol));
	SMP = smatrix_multiply(SMP, smatrix_layer(kh2,eps2,kx0));
	TMP = tmatrix_interface(eps2,eps1,kx0,pol);
	TMP = TMP * tmatrix_layer(kh1,eps1,kx0);
	TMP = TMP * tmatrix_interface(eps1,eps2,kx0,pol);
	TMP = TMP * tmatrix_layer(kh2,eps2,kx0);

		% matrices of the whole slab:
		% the first period contacts the surrounding medium:
	SM = smatrix_interface(eps,eps1,kx0,pol);
	SM = smatrix_multiply(SM, smatrix_layer(kh1,eps1,kx0));
	SM = smatrix_multiply(SM, smatrix_interface(eps1,eps2,kx0,pol));
	SM = smatrix_multiply(SM, smatrix_layer(kh2,eps2,kx0));
	TM = tmatrix_interface(eps,eps1,kx0,pol);
	TM = TM * tmatrix_layer(kh1,eps1,kx0);
	TM = TM * tmatrix_interface(eps1,eps2,kx0,pol);
	TM = TM * tmatrix_layer(kh2,eps2,kx0);
		% np-2 periods inside the slab
	for i = 2:(np-1)
		SM = smatrix_multiply(SM, SMP);
		TM = TM*TMP;
	end
		% the last period contacts the surrounding medium:
	SM = smatrix_multiply(SM, smatrix_interface(eps2,eps1,kx0,pol));
	SM = smatrix_multiply(SM, smatrix_layer(kh1,eps1,kx0));
	SM = smatrix_multiply(SM, smatrix_interface(eps1,eps2,kx0,pol));
	SM = smatrix_multiply(SM, smatrix_layer(kh2,eps2,kx0));
	SM = smatrix_multiply(SM, smatrix_interface(eps2,eps,kx0,pol));
	TM = TM * tmatrix_interface(eps2,eps1,kx0,pol);
	TM = TM * tmatrix_layer(kh1,eps1,kx0);
	TM = TM * tmatrix_interface(eps1,eps2,kx0,pol);
	TM = TM * tmatrix_layer(kh2,eps2,kx0);
	TM = TM * tmatrix_interface(eps2,eps,kx0,pol);

		% store power reflection coefficients
	reflection_ST(iwl,1) = abs(SM(2,2)^2);
	reflection_ST(iwl,2) = abs((TM(2,1)/TM(2,2))^2);
	%reflection_ST(iwl,3) = abs(SME(2,2)^2);
end

plot(wls, reflection_ST(:,1),'r');
hold on;
plot(wls, reflection_ST(:,2),'b');
hold off;
xlabel('\lambda');
ylabel('R');
ylim([0 1.1])

%}
%{
	% reflection from a homogeneous uniaxial slab
		% effectivie uniaxial medium slab
	SME = zeros(2,2);
	if strcmp(pol,'TE')
		kz1 = sqrt(epsR_o - kx0*kx0);
		kz2 = sqrt(eps - kx0*kx0);
	elseif strcmp(pol,'TM')
		kz1 = sqrt(eps_o - (eps_o/eps_e)*kx0*kx0)/eps_o;
		kz2 = sqrt(eps - kx0*kx0)/eps;
	end
	R = (kz2 - kz1)/(kz1 + kz2); % reflection coefficient at the boundary of the effective medium and the surrounding medium
	SME(1,1) = R*(1 - (1 - R^2)*exp(2*1i*kh*kz1)/(1 - R^2*exp(2*1i*kh*kz1)));
	SME(2,2) = SME(1,1);
	SME(1,2) = (1 - R^2)*exp(1i*kh*kz1)/(1 - R^2*exp(2*1i*kh*kz1));
	SME(2,1) = SME(1,2);
%}