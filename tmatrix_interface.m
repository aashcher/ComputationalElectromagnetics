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
% calculate 2-by-2 T-matrix of an interface between two homogeneous isotropic media
%% input:
% eps1, eps2: medium permittivities below and above the interface
% kx0: in-plane wavevector projection
% pol: polarization, either 'TE' or 'TM'
%% output:
% SMI: T-matrix of size 2-by-2
%% implementation:
function [TMI] = tmatrix_interface(eps1, eps2, kx0, pol)
    TMI = zeros(2,2);
    SMI = smatrix_interface(eps1, eps2, kx0, pol);
    tv = 1/SMI(2,1);
    TMI(2,2) = tv;
    TMI(2,1) = tv*SMI(2,2);
    TMI(1,2) = -tv*SMI(1,1);
    TMI(1,1) = SMI(1,2) - tv*SMI(2,2)*SMI(1,1);
end