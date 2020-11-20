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
% calculation of the band diangrams for 2D photonic crystals made of infinite
% circular rods with square or hexagonal lattice by the Fourier method
%% implementation:
%%
clc;

%% parameter definitions:

N = 15; % number of Fourier harmonics in one dimension
R = 0.4; % radius of cylinders (rods)
P = 1; % period (P > 2R)
kR = R/P; % radius-to-period ratio
eps1 = 14; % permittivity of cylinders
eps2 = 1; % permittivity of the surrounding medium
pol = 'TE'; % polarization
lattice_type = "square"; % "square" or "hexagonal"

kg = zeros(2,2); % matrix for the reciprocal lattice vectors

if strcmp(lattice_type,"square")
	kg(1,1) = 1; kg(1,2) = 0; % first reciprocal lattice vector
	kg(2,1) = 0; kg(2,2) = 1; % second reciprocal lattice vector
	f = pi*R^2/P^2; % volume fraction of cylinders
elseif strcmp(lattice_type,"hexagonal")
	kg(1,1) = 1; kg(1,2) = -1/sqrt(3);
	kg(2,1) = 0; kg(2,2) = 2/sqrt(3);
	f = (2*pi/sqrt(3))*R^2/P^2; % volume fraction of cylinders
end

%% test run:
%{
	% Bloch wavevector projections
kx0 = 0.3;
ky0 = 0.3;
	% calculate the reciprocal lattice vectors:
G = get_G(kx0,ky0,kg,N,N);
	% calculate the permittivity Fourier matrix:
FE = get_eps_Fourier_cyl(N,N,f,kR,kg,eps1,eps2);
	% solve the eigenvalue problem:
M = get_2Deig_matrix(N,N,G,FE,pol);
W = sort(sqrt(eig(M))); % vector of eigen frequencies
disp(W);

return;
%}
%% band diagram calculation and visualization:
	% set the number of bands to store:
n_band = 5;
	% initialization of the band diagram calculation:
  % set numbers of points to calculate and coordinates of corner points
  % in the reciprocal space
if strcmp(lattice_type,"square")
	n12 = 30; x1 = 0; y1 = 0; % Gamma point
	n23 = 30; x2 = 0.5; y2 = 0.5; % M point
	n31 = 30; x3 = 0.5; y3 = 0; % X point
elseif strcmp(lattice_type,"hexagonal")
	n12 = 30; x1 = 0; y1 = 0; % Gamma point
	n23 = 30; x2 = 0.5; y2 = 0.5/sqrt(3); % M point
	n31 = 30; x3 = 2/3; y3 = 0; % K point
end
% {
%% band diagram calculation
nn = n12 + n23 + n31;
FE = get_eps_Fourier_cyl(N,N,f,kR,kg,eps1,eps2);

dataDD = zeros(n_band+1, nn); % matrix to store the data

  % line between points 1 and 2
for i=1:n12
    kx0 = (i-0.5)/n12 * (x2-x1) + x1;
    ky0 = (i-0.5)/n12 * (y2-y1) + y1;
    G = get_G(kx0,ky0,kg,N,N);
    M = get_2Deig_matrix(N,N,G,FE,pol);
    W = sort(sqrt(eig(M)));
    dataDD(1,i) = sqrt((kx0-x1)^2 + (ky0-y1)^2);
    for k=1:5
        dataDD(k+1,i) = W(k);
    end
end
kk = sqrt((x2-x1)^2 + (y2-y1)^2);
  % line between points 2 and 3
for i=1:n23
    kx0 = ((i-0.5)/n23) * (x3-x2) + x2;
    ky0 = ((i-0.5)/n23) * (y3-y2) + y2;
    G = get_G(kx0,ky0,kg,N,N);
    M = get_2Deig_matrix(N,N,G,FE,pol);
    W = sort(sqrt(eig(M)));
    dataDD(1,n12+i) = kk + sqrt((kx0-x2)^2+(ky0-y2)^2);
    for k=1:5
        dataDD(k+1,n12+i) = W(k);
    end
end
kk = kk + sqrt((x3-x2)^2 + (y3-y2)^2);
  % line between points 3 and 1
for i=1:n31
    kx0 = ((i-0.5)/n31) * (x1-x3) + x3;
    ky0 = ((i-0.5)/n31) * (y1-y3) + y3;
    G = get_G(kx0,ky0,kg,N,N);
    M = get_2Deig_matrix(N,N,G,FE,pol);
    W = sort(sqrt(eig(M)));
    dataDD(1,n12+n23+i) = kk + sqrt((kx0-x3)^2 + (ky0-y3)^2);
    for k=1:5
        dataDD(k+1,n12+n23+i) = W(k);
    end
end

% plot the band diagram:
for k = 1:n_band
    plot(dataDD(1,:),dataDD(1+k,:));
    hold on;
end
hold off;
xlabel('k-vector trajectory')
ylabel('\omega\Lambda/2\pi c')
xticks([0 (sqrt((x2-x1)^2+(y2-y1)^2)) kk (kk+sqrt((x1-x3)^2 + (y1-y3)^2))])
xticklabels({'\Gamma','M','X/K','\Gamma'})
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%

%% calculate the matrix of reciprocal lattice vectors
% INPUT:
% kx0, ky0: Bloch wavevector projections (dimensionless, normalized by the period)
% kg: reciprocal lattice vectors (2x2 matrix)
% N1, N2: - numbers of Fourier harmonics in each dimension
function [G] = get_G(kx0, ky0, kg, N1, N2)
	G = zeros(3,N1*N2);
	for i1 = 1:N1
		for i2 = 1:N2
			G(1,(i1-1)*N2+i2) = kx0 + (i1-1-floor(N1/2))*kg(1,1) + (i2-1-floor(N2/2))*kg(2,1);
			G(2,(i1-1)*N2+i2) = ky0 + (i1-1-floor(N1/2))*kg(1,2) + (i2-1-floor(N2/2))*kg(2,2);
			G(3,(i1-1)*N2+i2) = sqrt((G(1,(i1-1)*N2+i2))^2 + (G(2,(i1-1)*N2+i2))^2);
		end
	end
end

%% calculate matrix of the Fourier components of the inverse permittivity
% INPUT:
% N1, N2: number of Fourier harmonics in each direction
% f: volume filling factor (rod cross section area divided by the period area)
% kR: rod radius divided by period
% kg: reciprocal lattice vectors (2x2 matrix)
% eps1,2: permittivities of rods (cylinders) and the surrounding medium
function [FE] = get_eps_Fourier_cyl(N1, N2, f, kR, kg, eps1, eps2)
	FE = zeros(2*N1,2*N2);
	tc = f*(1/eps1 - 1/eps2);
	for i1 = 1:(2*N1)
		m1 = i1 - N1;
		for i2 = 1:(2*N2)
			m2 = i2 - N2;
			if (m1 == 0) && (m2 == 0)
				FE(i1,i2) = tc + 1/eps2;
			else
				RG = 2*pi*kR*sqrt((kg(1,1)*m1 + kg(2,1)*m2)^2 + (kg(1,2)*m1 + kg(2,2)*m2)^2);
				FE(i1,i2) = 2*tc*besselj(1,RG)/RG;
			end
		end
	end
end

%% calculate matrix of the Fourier components of the 
% N1, N2 - neumber of Fourier harmonics in each direction
% kG - matrix of reciprocal lattice vectors (including Bloch vawevector)
% FE - matrix of the Fourier components of the inverse permittivity
% pol - polarization, TE or TM
function [M] = get_2Deig_matrix(N1, N2, kG, FE, pol)
	M = zeros(N1*N2,N1*N2);
	if strcmp(pol,'TE')
		for i1 = 1:N1
			for i2 = 1:N2
				i = (i1-1)*N2 + i2;
				for ii1 = 1:N1
					for ii2 = 1:N2
						ii = (ii1-1)*N2 + ii2;
						M(i,ii) = kG(3,i)*kG(3,ii)*FE(ii1-i1+N1,ii2-i2+N2);
					end
				end
			end
		end
	else
		for i1 = 1:N1
			for i2 = 1:N2
				i = (i1-1)*N2 + i2;
				for ii1 = 1:N1
					for ii2 = 1:N2
						ii = (ii1-1)*N2 + ii2;
						M(i,ii) = (kG(1,i)*kG(1,ii) + kG(2,i)*kG(2,ii))*FE(ii1-i1+N1,ii2-i2+N2);
					end
				end
			end
		end
	end
end

% END %