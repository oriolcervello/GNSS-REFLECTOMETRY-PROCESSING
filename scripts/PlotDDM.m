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
% FOR DDM
% -----------------------
%%
%OPEN FILE
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\TEST_beam_1_s_4_d/PeaksIteration1380.bin','r');
% fileID = fopen('prn_L1CA_32.bin','r');
%%
%READING COMPLEX SIGNALS
A = fread(fileID,[2 1000000000000000],'float32');  % Per senyals complexes
A = A(1,:)+1j*A(2,:);
fclose(fileID);

%%
samplesToSave=311;
numOfFFT=1000;
ddmQuant=21;

B = reshape(A,[samplesToSave ddmQuant*numOfFFT]);
% B = reshape(A,[ddmQuant samplesToSave*numOfFFT]);

%%

% C = B(:,1:1000:end);
C = B(1:311,1:1000:end);


%%
Span = 11;
Resolucio = 100; % Hz
F_doppler = 517; % Hz
F_ = F_doppler-(Span-1)/2*Resolucio:Resolucio:F_doppler+(Span-1)/2*Resolucio;

%%
 ddmToView=1;
% 
C=B(:,ddmToView:numOfFFT:ddmQuant*numOfFFT);

C=abs(C').^2;
% C=abs(C);
 figure
% mesh(F_,1:samplesToSave,C);
mesh(C)


%%

D = abs(B(:,1:numOfFFT:ddmQuant*numOfFFT)).^2;
for i=2:1000
    D = D+abs(B(:,i:numOfFFT:ddmQuant*numOfFFT)).^2;
end

 figure
mesh(D)


%%
figure
for i=1:numOfFFT
    C=B(:,i:numOfFFT:ddmQuant*numOfFFT);
    D = C;
    C=abs(C).^2;
    
    [a,b] = max(max(C));
    [~,c] = max(max(C'));
    
%     mesh(F_,1:samplesToSave,C);
    mesh(C)
    
    valor(i) = a;
    delay(i) = c;
    dopp(i) = b;
    
    valorc(i) =D(c,b);
    
    pause(0.1);
end
valor = valor/32768;
valorc = valorc/32768;


%% Evolucio temporal


Fs = 1/1e-3;
temps = linspace(0,length(valorc)/Fs,length(valorc));

figure, 
        scatter(real(valorc),imag(valorc),100,temps,'filled')
        grid minor
        AX = colorbar;
        AX.Label.String = 'Temps [s]';
        xlabel('I Ampltiud [A.U.]')
        ylabel('Q Ampltiud [A.U.]')
figure, 
        hold on;
        plot(temps,real(valorc),'b')
        plot(temps,imag(valorc),'r')
        legend('I','Q')
        grid minor
        xlabel('Temps [s]')
        ylabel('Ampltiud [A.U.]')


figure, plot(temps,10*log10(valor))
        grid minor
        xlabel('Temps [s]')
        ylabel('Power [dB A.U.]')

figure, stairs(temps,delay)
        grid minor
        xlabel('Temps [s]')
        ylabel('Samples')
figure, stairs(temps,F_(dopp)/1000)
        grid minor
        xlabel('Temps [s]')
        ylabel('Doppler [kHz]')
        
%% Histogrames 

figure; hist(valor)
figure; hist(delay)
figure; hist(F_(dopp))

%% Covariances

aux = cov([real(valorc); fliplr(real(valorc))],1);
figure, imagesc(temps,temps,(aux));
        set(gca,'YDir','normal');
        xlabel('Time')
        ylabel('Time')
        colorbar
        
aux = cov([imag(valorc); fliplr(imag(valorc))],1);
figure, imagesc(temps,temps,(aux));
        set(gca,'YDir','normal');
        xlabel('Time')
        ylabel('Time')
        colorbar        

aux = cov([real(valorc); imag(valorc)],1);
figure, imagesc(temps,temps,(aux));
        set(gca,'YDir','normal');
        xlabel('Time')
        ylabel('Time')
        colorbar

aux = cov([valorc;fliplr(valorc)],1);
figure, imagesc(temps,temps,real(aux));
        set(gca,'YDir','normal');
        xlabel('Time')
        ylabel('Time')
        colorbar
        
figure, imagesc(temps,temps,imag(aux));
        set(gca,'YDir','normal');
        xlabel('Time')
        ylabel('Time')
        colorbar        
       



%%
fclose(fileID);
clear all;
