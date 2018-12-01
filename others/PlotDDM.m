% -----------------------
% FOR DDM
% -----------------------
%%
%OPEN FILE
fileID = fopen('C:\Users\ori\Documents\Estudis\uni\SA\tfg\CUDA_repos\PROJECT\cufft2\results\PeaksIteration0.bin','r');
% fileID = fopen('prn_L1CA_32.bin','r');
%%
%READING COMPLEX SIGNALS
A = fread(fileID,[2 1000000000000000],'float32');  % Per senyals complexes
A = A(1,:)+1j*A(2,:);

%%
samplesToSave=311;
numOfFFT=80;
ddmQuant=51;

B = reshape(A,[samplesToSave ddmQuant*numOfFFT]);

%%
Span = 51;
Resolucio = 25; % Hz
F_doppler = 517; % Hz
F_ = F_doppler-(Span-1)/2*Resolucio:Resolucio:F_doppler+(Span-1)/2*Resolucio;

%%
% ddmToView=17;
% 
% C=B(:,ddmToView:numOfFFT:ddmQuant*numOfFFT);
% C=abs(C).^2;
% % C=abs(C);
% %  figure
% % mesh(F_,1:samplesToSave,C);

%%

for i=1:numOfFFT
    C=B(:,i:numOfFFT:ddmQuant*numOfFFT);
    D = C;
    C=abs(C).^2;
    
    [a,b] = max(max(C));
    [~,c] = max(max(C'));
    
    mesh(F_,1:samplesToSave,C);
    
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
