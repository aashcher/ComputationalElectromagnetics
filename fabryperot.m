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
% calculate a reflection map of a Fabry-Perot resonator
%% implementation:

h = 0.5;
eps_s = 1; % permittivity of the surrounding medium
eps_l = 4; %
theta = pi/10;

kx = sqrt(eps_s)*sin(theta);
kz_s = sqrt(eps_s - kx^2);
kz_l = sqrt(eps_l - kx^2);
Re = ((kz_l - kz_s)/(kz_l + kz_s))^2;
Rh = ((eps_s*kz_l - eps_l*kz_s)/(eps_s*kz_l + eps_l*kz_s))^2;

wl = linspace(0.4,0.8,100);
Te = (1-Re)^2 * ((1-Re)^2 + 4*Re*sin(2*pi*h./wl*kz_l).^2).^(-1);
Th = (1-Rh)^2 * ((1-Rh)^2 + 4*Rh*sin(2*pi*h./wl*kz_l).^2).^(-1);
fr = 1./wl;

%plot(fr,Te);

%% reflection map
% {
[wl,theta] = meshgrid(0.4:0.002:0.8,0:0.5:89);
kx = sin((pi/180)*theta);
kz_s = sqrt(eps_s - kx.^2);
kz_l = sqrt(eps_l - kx.^2);
Re = ((kz_l - kz_s)./(kz_l + kz_s)).^2;
%Rh = ((eps_s*kz_l - eps_l*kz_s)./(eps_s*kz_l + eps_l*kz_s)).^2;
Te = (1-Re).^2 .* ((1-Re).^2 + 4*Re.*sin((2*pi*h)./wl.*kz_l).^2).^(-1);
%Th = (1-Rh)^2 * ((1-Rh)^2 + 4*Rh*sin(2*pi*h./wl*kz_l).^2).^(-1);

surf(theta,wl,Te,'EdgeColor','none');
view(2);
xlabel('\theta');
ylabel('\lambda');
