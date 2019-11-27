%==========================================================================
% Author: Oriol Cervelló (oriol.cn [at] protonmail.com) 
%==========================================================================
% License: GNU GPLv3.0
% Copyright (C) 2019  Oriol Cervelló
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%==========================================================================
% -----------------------
% FOR READING MAX FILE
% -----------------------
%%
%OPEN FILE
% fileID = fopen('Maximums.bin','r');
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\TEST_beam_1_s_2_d/Maximums300.bin','r');
% fileID = fopen('prn_L1CA_32.bin','r');
%%
%READING FILE
A = fread(fileID,[5 1000000000000000],'single');  % Per senyals complexes

%%
%WHAT IS EACH COLUM
posOfMax=A(1,:);
maxValue=A(2,:);
meanValue=A(3,:);
stdValue=A(4,:);
dopplerFreq=A(5,:);
%%
% FIND MAX FOR 1 FFT
numOfFFT=1000;
averageIncoherent=10;

FFTtoView=1;
aux=A(:,FFTtoView:numOfFFT/averageIncoherent:length(maxValue));

[value,pos]=max(aux(2,:));
valueCenteredForSavingPeak=aux(:,pos);
%%
fclose(fileID);
clear all;

%%
for i=1:100

    aux = A(:,i:1000/10:end);
    [~,pos] = max(aux(2,:));
%     aux(1,pos)
    caca(i) = aux(2,pos)
end

figure, hist(caca.^2)
