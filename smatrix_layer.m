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
% calculate 2-by-2 S-matrix of a homogeneous layer
%% input:
% kh: layer thickness multiplied by the vacuum wavenumber
% kx0: in-plane wavevector projection 
% eps: layer permittivity
%% output:
% SML: S-matrix of size 2-by-2
%% implementation:
function [SML] = smatrix_layer(kh, eps, kx0)
    SML = zeros(2,2);
    tv = exp(1i*kh*sqrt(eps - kx0^2));
    SML(2,1) = tv;
    SML(1,2) = tv;
end