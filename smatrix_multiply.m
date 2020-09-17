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
% S-matrix multiplication
%% input:
% SM1, SM2: low and upper S-matrices of size 2-by-2
%% output:
% SM: composition S-matrix of size 2-by-2
%% implementation:
function [SM] = smatrix_multiply(SM1, SM2)
    tv = 1/(1 - SM1(2,2)*SM2(1,1));
    SM(2,1) = tv*SM1(2,1)*SM2(2,1);
    SM(1,1) = SM1(1,1) + tv*SM1(2,1)*SM1(1,2)*SM2(1,1);
    SM(1,2) = tv*SM2(1,2)*SM1(1,2);
    SM(2,2) = SM2(2,2) + tv*SM2(1,2)*SM2(2,1)*SM1(2,2);
end