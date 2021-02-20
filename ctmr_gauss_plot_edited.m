function [c_h] = ctmr_gauss_plot_edited(cortex,elecmatrix,weights,cax,addl,CM,gsp)
% function [c_h]=ctmr_gauss_plot(cortex,elecmatrix,weights)
%
% projects electrode locationsm (elecmatrix) onto their cortical spots in 
% the left hemisphere and plots about them using a gaussian kernel
% for only cortex use:
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)

%     Copyright (C) 2009  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   Version 1.1.0, released 26-11-2009
%
% Modified in 2016 by Liberty Hamilton and then by Jon Kleen in 2017-2021
% for omni-planar and surface casting of epileptiform activity (OPSCEA).


if isempty(elecmatrix); 
    elecmatrix = [0 0 0];
end
if isempty(weights); 
    weights = zeros(size(elecmatrix,1),1);
end

if ~exist('addl','var'); addl=0; end

if exist('CM','var') && ~isempty(CM); cm=CM; else cm=cmSz3D; end

brain=cortex.vert;

if length(weights)~=length(elecmatrix(:,1))
    error('You sent a different number of weights than electrodes in elecmatrix (perhaps a whole matrix instead of vector)')
end

if ~exist('gsp','var') || isempty(gsp); gsp=10; end %default 10

c=zeros(length(cortex(:,1)),1);
for i=1:length(elecmatrix(:,1))
    b_z=abs(brain(:,3)-elecmatrix(i,3));
    b_y=abs(brain(:,2)-elecmatrix(i,2));
    b_x=abs(brain(:,1)-elecmatrix(i,1));
    d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian
    c=c+d';
end

c_h=tripatch(cortex, 'nofigure', c');
if ~addl; shading interp; end
a=get(gca);

d=a.CLim;
if exist('cax','var') && ~isempty(cax); set(gca,'CLim',[cax(1) cax(2)]); 
else set(gca,'CLim',[-max(abs(d)) max(abs(d))]); 
end
colormap(gca,cm)
material dull;
axis off
litebrain('a',.1) %just so it is visible; will replace lighting soon