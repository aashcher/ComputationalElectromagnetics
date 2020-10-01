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
% plot the dispersion function for a two dielectric layer 1D photonic crystal
% and calculate mode propagation constants
%% implementation:
%%
clc;
%% parameters of the photonic crystal:
eps1 = 1; % layer 1 permittivity
eps2 = 8; % layer 2 permittivity
d1 = 0.4; % layer 1 thickness
d2 = 0.5; % layer 2 thickness
L = d1 + d2; % period
wl = 0.5; % wavelength
kd1 = 2*pi*d1/wl; % normalized layer 1 thickness
kd2 = 2*pi*d2/wl; % normalized layer 2 thickness

kb = 0.1; % normalized Bloch wavenumber
ck = cos(2*pi*kb); % right-hand constant in the dispersion equation

%% plot the dispersion function:
plot_dispersion_func(eps1,eps2,kd1,kd2,ck,'TE');

%% roots of the dispersion equation:

	% dispersion as a function of the squared propagation constant:
fun = @(x) dispersion_PhC1D(x,eps1,eps2,kd1,kd2,ck,'TE');

	% roots in the region [eps1,eps2]:
nr = 2*(floor(kd2*sqrt(eps2-eps1)/pi)+1); % number of initial points
betas = []; % vector of squared propagation constants
for i = 1:nr
	b_ini = eps2 - (eps2-eps1)*((i-0.5)/nr).^2; % initial point for fzero()
	tb = fzero(fun,b_ini);
	if numel(find(abs(betas-tb) < 1e-3)) == 0 && tb > eps1
		betas = [betas, tb]; % if a new root is found, then store it
	end
end
betas = sort(betas,'descend');

	% roots in the region < eps1:
for k = 0:9
		% estimate the next region containing one or two roots:
	b1 = min(eps1,eps2) - (k*pi/max(kd1,kd2))^2;
	b2 = min(eps1,eps2) - ((k+1)*pi/max(kd1,kd2))^2;
	for kk = 1:4 % take four initial point in this region
		b_ini = b1 - ((kk-0.5)/4)*(b1-b2);
		tb = fzero(fun,b_ini);
		if tb <= b1 && tb >= b2 && abs(tb - betas(numel(betas))) > 1e-10 && abs(tb - betas(numel(betas)-1)) > 1e-10
			betas = [betas, tb]; % if a new root is found, then store it
		end
	end
end
disp(betas);

%% dispersion function:
function f = dispersion_PhC1D(beta2, eps1, eps2, kd1, kd2, C, pol)
	k1d = kd1*sqrt(eps1 - beta2);%^2
	k2d = kd2*sqrt(eps2 - beta2);%^2
	if strcmp(pol,'TE')
		tvar = k1d./k2d;
	elseif strcmp(pol,'TM')
		tvar = (eps2/eps1)*k1d./k2d;
	end
	f = cos(k1d).*cos(k2d) - 0.5*(tvar + 1./tvar).*sin(k1d).*sin(k2d) - C;
end

%% plot the dispersion function:
function plot_dispersion_func(eps1, eps2, kd1, kd2, C, pol)
	bmin = min(eps1,eps2);
	bmax = max(eps1,eps2);
	x_min = -3*bmax;
	x_max = 1.1*bmax;
	beta2 = x_min : 0.001*(x_max-x_min) : x_max;
	func = dispersion_PhC1D(beta2, eps1, eps2, kd1, kd2, 0, pol);
	
	plot(beta2,func);

	line([-5*bmax,1.1*bmax],[1,1],'Color','black','LineStyle','--');
	line([-5*bmax,1.1*bmax],[-1,-1],'Color','black','LineStyle','--');
	line([-5*bmax,1.1*bmax],[C,C],'Color','red','LineStyle','-');

	line([bmin,bmin],[-4.5,4.5],'Color','black','LineStyle','--');
	line([bmax,bmax],[-4.5,4.5],'Color','black','LineStyle','--');

%	for k = 0:25
%		bb = (eps1) - (k*pi/kd2)^2;
%		line([bb,bb],[-4.5,4.5],'Color','green','LineStyle','--');
%	end
	
	xlim([x_min x_max]);
	ylim([-4.5 4.5]);

	xlabel('\beta^2');
	ylabel('f');
end


%%% END OF FILE