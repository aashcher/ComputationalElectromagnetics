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
% calculate 2-by-2 S-matrix of an interface between two homogeneous isotropic media
%% input:
% eps1, eps2: medium permittivities below and above the interface
% kx0: in-plane wavevector projection
% pol: polarization, either 'TE' or 'TM'
%% output:
% SMI: S-matrix of size 2-by-2
%% implementation:
function [SMI] = smatrix_interface(eps1, eps2, kx0, pol)
    SMI = zeros(2,2);
    if strcmp(pol,'TE')
        kz1 = sqrt(eps1 - kx0^2);
        kz2 = sqrt(eps2 - kx0^2);
    elseif strcmp(pol,'TM')
        kz1 = sqrt(eps1 - kx0^2)/eps1;
        kz2 = sqrt(eps2 - kx0^2)/eps2;
    end
    SMI(1,1) = (kz1 - kz2)/(kz1 + kz2);
    SMI(1,2) = 1 - SMI(1,1);
    SMI(2,2) = -SMI(1,1);
    SMI(2,1) = 1 + SMI(1,1);
end


