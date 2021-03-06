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
% compute and plot the eigen modes for a two dielectric layer 1D photonic crystal
%% implementation:
clc;
%%
d1 = 0.4; % first layer thickness
d2 = 0.6; % second layer thickness
eps1 = 1; % first layer permittivity
eps2 = 5; % second layer permittivity
L = d1 + d2; % period
wl = 0.5; % wavelength
kd1 = 2*pi*d1/wl; % normalized layer 1 thickness
kd2 = 2*pi*d2/wl; % normalized layer 2 thickness
kL = 2*pi*L/wl; % normalized period

kb = 0; % normalized Bloch wavenumber (sin(theta_inc))
ck = cos(kb*kL); % right-hand side constant in the dispersion equation

%% calculate propagation constants
	% dispersion as a function of the squared propagation constant:
fun = @(x) dispersion_PhC1D(x,eps1,eps2,kd1,kd2,ck,'TE');
	% find roots in the region [eps1,eps2]:
nr = 2*(floor(kd2*sqrt(eps2-eps1)/pi)+1); % number of initial points
betas = []; % vector of squared propagation constants
options = optimset('TolX',1e-14);
for i = 1:nr
	b_ini = eps2 - (eps2-eps1)*((i-0.5)/nr).^2; % initial point for fzero()
	tb = fzero(fun,b_ini,options);
	if numel(find(abs(betas-tb) < 1e-3)) == 0 && tb > eps1
		betas = [betas, tb]; % if a new root is found, then store it
	end
end
betas = sort(betas,'descend');
	% find roots in the region < eps1:
for k = 0:9
		% estimate the next region containing one or two roots:
	b1 = min(eps1,eps2) - (k*pi/max(kd1,kd2))^2;
	b2 = min(eps1,eps2) - ((k+1)*pi/max(kd1,kd2))^2;
	for kk = 1:4 % take four initial point in this region
		b_ini = b1 - ((kk-0.5)/4)*(b1-b2);
		tb = fzero(fun,b_ini);
		if tb <= b1 && tb >= b2 && abs(tb - betas(numel(betas))) > 1e-10
			betas = [betas, tb]; % if a new root is found, then store it
		end
	end
end

disp(betas);
betas = sqrt(betas);

%% calculate and plot fields
m_ind = 2; % number of mode to plot
kappa1 = get_kappa(betas(m_ind),eps1);
kappa2 = get_kappa(betas(m_ind),eps2);
	% z points in a single period
z = 0 : 0.001*kL : kL;
	% calculate the modal field amplitudes for the set of points
[F,G] = get_fmodal(z,kb,kappa1,kappa2,kd1,kd2);
	% plot three periods for real and imaginary parts of the modal function
z = z*(1/kL) - 0.5*(d1/L);
plot(z,real(F),'Color','r');
hold on;
plot(z+1,real(F),'Color','r');
hold on;
plot(z-1,real(F),'Color','r');
hold on;
plot(z,imag(F),'Color','b');
hold on;
plot(z+1,imag(F),'Color','b');
hold on;
plot(z-1,imag(F),'Color','b');
hold off;

return;

%% calculate kappa
function kappa = get_kappa(beta, eps)
	kappa = sqrt(eps - beta.^2);
	ind = angle(kappa) < -1e-10;
	kappa(ind) = -kappa(ind);
end

%% calculate modal field
	% F - E_y field amplitude (for the TE polarization)
	% G = dF/dz
function [F,G] = get_fmodal(z, k0, kappa1, kappa2, d1, d2)
	ind1 = z <= d1;
	ind2 = ~ind1;
	F = z;
	G = z;

	tau = exp(1i*k0*(d1 + d2));
	ed1 = exp(1i*kappa1*d1);
	ed2 = exp(1i*kappa2*d2);
	kk = kappa1/kappa2;
	
	a1p = 1;
	a1m = (1+kk)*(tau*ed1*ed2-1)/(1-kk)/(1-tau*ed2/ed1);
	a2p = 0.5*( a1p*(1+kk)*ed1 + a1m*(1-kk)/ed1 );
	a2m = 0.5*( a1p*(1-kk)*ed1 + a1m*(1+kk)/ed1 );

	F(ind1) = ( a1p * exp(1i*kappa1*z(ind1)) + a1m * exp(-1i*kappa1*z(ind1)) );
	F(ind2) = ( a2p * exp(1i*kappa2*(z(ind2)-d1)) + a2m * exp(-1i*kappa2*(z(ind2)-d1)) );
	G(ind1) = 1i*( (k0+kappa1) * a1p * exp(1i*kappa1*z(ind1)) ...
										+ (k0-kappa1) * a1m * exp(-1i*kappa1*z(ind1)) );
	G(ind2) = 1i*( (k0+kappa2) * a2p * exp(1i*kappa2*(z(ind2)-d1)) ...
										+ (k0-kappa2) * a2m * exp(-1i*kappa2*(z(ind2)-d1)) );
end

%% dispersion function:
function f = dispersion_PhC1D(beta2, eps1, eps2, kd1, kd2, C, pol)
	kp1 = sqrt(eps1 - beta2);
	kp2 = sqrt(eps2 - beta2);
	k1d = kd1*kp1;
	k2d = kd2*kp2;
	if strcmp(pol,'TE')
		tvar = kp1./kp2;
	elseif strcmp(pol,'TM')
		tvar = (eps2/eps1)*kp1./kp2;
	end
	f = cos(k1d).*cos(k2d) - 0.5*(tvar + 1./tvar).*sin(k1d).*sin(k2d) - C;
end

%% END