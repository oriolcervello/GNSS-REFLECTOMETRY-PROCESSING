%==========================================================================
% Author: Oriol Cervelló (oriol.cn [at] protonmail.com) 
%==========================================================================
% License: GPL-3.0-only
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
% figure
% clear all;
% numfitxers= 574;
carpeta='YANCO_DRY_beam_4_s_2_d';
numOfFFT=1000;
% averageIncoherent=250;
% averageIncoherent=100;
averageIncoherent=10;
% averageIncoherent=1;

% FFT_size = 32768;
FFT_size = 32736;

a=dir(['C:\GNSS-REFLECTOMETRY-PROCESSING\results\YANCO\' num2str(FFT_size) '\incoh_' num2str(averageIncoherent) '\' carpeta ]);
numfitxers =size(a,1)-2-1


cont=1;
maximums=zeros(5,numfitxers*numOfFFT/averageIncoherent);
[a,b] = size(maximums);

for i=1:numfitxers
nomfitxer=['C:\GNSS-REFLECTOMETRY-PROCESSING\results\YANCO\'  num2str(FFT_size) '\incoh_' num2str(averageIncoherent) '\' carpeta '\Maximums' num2str(i-1) '.bin']
fileID = fopen(nomfitxer,'r');
A = fread(fileID,[5 1000000000000000],'single');

if i==1
    aux = length(A);
    if aux == 250
        break;
    end
end

% DDMs
% % % % posOfMax=A(1,:);
% % % maxValue=A(2,:);
% % % % meanValue=A(3,:);
% % % % stdValue=A(4,:);
% % % % dopplerFreq=A(5,:);
% % % 
% % % 
% % % for FFTtoView=1:numOfFFT/averageIncoherent
% % % aux=A(:,FFTtoView:numOfFFT/averageIncoherent:length(maxValue));
% % % 
% % % % [value,pos]=max(aux(2,:));
% % % pos = 11;
% % % valueCenteredForSavingPeak=aux(:,pos);
% % % 
% % % maximums(:,cont)=valueCenteredForSavingPeak(:);
% % % 
% % % 
% % % cont=cont+1;
% % % end


% Waveform
maximums(:,cont:cont-1+numOfFFT/averageIncoherent)=A(:,:);
cont=cont+numOfFFT/averageIncoherent;



fclose(fileID);

end
%%

% figure; 
hold on;
    %scatter(linspace(0,b*averageIncoherent*1e-3,b),maximums(1,:)*1/(32*1.023e6*.1e-3/10230),100,maximums(2,:),'filled')
%     scatter(linspace(0,b*averageIncoherent*1e-3,b),maximums(1,:)*1/(32*1.023e6*.1e-3/10230),5,'filled')
%     scatter(linspace(0,b*averageIncoherent*1e-3,b),maximums(2,:),5,'filled')

% scatter(linspace(ref,ref+b*averageIncoherent*1e-3,b),maximums(1,:)*1/(32*1.023e6*.1e-3/10230),5,'filled')
scatter(linspace(ref,ref+b*averageIncoherent*1e-3,b),maximums(1,:)*1/(32*1.023e6*.1e-3/10230),20,maximums(2,:),'filled')
ref =ref+b*averageIncoherent*1e-3;
    
%%
% figure
%  scatter(linspace(0,b*10e-3,b),(DIRECTA(1,:)-REFLECTIDA(1,:))*1/(32*1.023e6/3e8),100,REFLECTIDA(2,:),'filled')
