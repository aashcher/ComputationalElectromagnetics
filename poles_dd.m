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
% functions based on the Newton divided differences to calculate comples
% poles of meromorphic functions, and examples
% for a reference, see
% [1] A.V. Tishchenko, M. Hamdoun, Opt. Quant. Electron. 35, 475 (2003)
%% implementation:
%%
clc;

n_poles = 2; % number of poles to search for
max_pow = 5; % maximum power of the Taylor series for the function f(z)P(z)

fun = @(x) test_function(x); % input function
	% calculate poles and pole amplitudes
[zp,A] = poles(fun,[-6.5,1.5],n_poles,max_pow);
for i = 1:n_poles
	fprintf("pole %d: z = (%e,%e), A = (%e,%e)\n", ...
		i,real(zp(i)), imag(zp(i)), real(A(i)), imag(A(i)));
end

%% functions %%
	% test meromorphic function
function f = test_function(z)
	f = 5 + (z+2).*(z+3)./(z-1+(0e-8)*1i)./(z+5);
end

	% pole calculation function:
function [zp,A] = poles(fun, z, n_poles, max_pow)
		% INPUT:
		% fun: reference to a funciton under consideration
		% z: 2-by-1 array of with real bounds of a search interval
		% n_poles: integer number of poles to search for
		% max_pow: maximum power of the Taylor series for the function f(z)P(z)
		% OUTPUT:
		% zp: vector of calulated n_poles poles
		% A: vector of calulated n_poles pole amplitudes
	if n_poles == 1
		zz = z(1) + 1e-10 + ((z(2)-z(1))/(max_pow))*(0 : max_pow);
		fp = fun(zz);

		Z = zz - zz';
		Z(1:max_pow+2:end) = 1;
		Pz = prod(Z,1);
		zp = sum(fp.*zz./Pz) / sum(fp./Pz);

		A = -sum(fp./Pz) * prod(zp-zz);
	else
		Nz = n_poles + max_pow;
		dz = (z(2)-z(1))/(Nz-1);
		nz = 1:1:n_poles;

		M = zeros(n_poles, n_poles);
		b = zeros(n_poles, 1);
			% calculte matrix elements
		for i = 1:n_poles
			zz = z(1) + 1e-5 + (i-1)*dz + dz*(0 : max_pow);
			fp = fun(zz);
			Z = zz - zz';
			Z(1:max_pow+2:end) = 1;
			Pz = prod(Z,1);
			b(i) = -sum(fp.*(zz.^(n_poles))./Pz);
			M(i,:) = (sum(fp.*(zz.^((0:n_poles-1)'))./Pz,2))';
		end
			% find coefficients of the polynomial
		p = M\b;
			% fill the companion matrix
		C = zeros(n_poles, n_poles);
		C(2:n_poles+1:end) = 1;
		C(:,n_poles) = -p;
			% find poles
		zp = eig(C);
			% calculate amplitudes:
		A = zeros(n_poles,1);
		Zp = zp' - zp;
		Zp(1:n_poles+1:end) = 1;
		PPz = prod(Zp,1);
		for k = 1:n_poles
			A(k) = (prod(zp(k)-zz)/PPz(k)) * sum(fp./(zp(k)-zz).*prod(zz-zp,1)./Pz);
		end
	end
end

% END %