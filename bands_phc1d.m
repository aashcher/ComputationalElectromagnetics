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
% calculate and plot the dispersion diagram for a two layer 1D photonic
% crystal in case of the normal to layer propagation
%% implementation:
	% photonic crystal parameters:
eps1 = 1; % permittivity of the first layer
eps2 = 4; % permittivity of the second layer
period = 0.2; % period in length units
d1 = 0.5*period; % thickness of the first layer in length units 
k = 0.1; % Bloch wavevector
	% numerical parameters
n = 11; % number of positive Fourier orders
N = 2*n; % number of positive Fourier orders for epsilon decomposition

nrow = linspace(-n,n,2*n+1);
Nrow = linspace(-N,N,2*N+1);
Nrow(N+1) = 1;

	% Fourier vector of the permittivity:
fv_eps = (eps2 - eps1) * sin((pi*d1/period)*Nrow) ./ ((pi)*Nrow);
fv_eps(N+1) = (eps1*d1 + eps2*(period - d1)) / period;
	% epsilon Fourier matrix
FM = toeplitz(fv_eps(N+1:2*N+1), fv_eps(N+1:-1:1));

nb = 4; % number of bands to calculate
kv = 0:0.001:0.5; % Bloch k-vectors
bands = zeros(numel(kv),nb);

	% loop over k-space points:
for ik = 1:numel(kv)
	A = diag((kv(ik) + nrow).^2);
	[V,D] = eig(A,FM);
	d = sort(real(sqrt(diag(D))),'ascend');
	bands(ik,:) = d(1:nb);
end
	% plot band diagram
for ib = 1:nb
	plot(kv, bands(:,ib), 'b');
	hold on;
end
hold off;
xlabel('k_x\Lambda/2\pi');
ylabel('\Lambda/\lambda');

return;

%% end %%