% -----------------------
% FOR READING MAX FILE
% -----------------------
%%
%OPEN FILE
fileID = fopen('Maximums.bin','r');
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
numOfFFT=100;
averageIncoherent=10;

FFTtoView=1;
aux=A(:,FFTtoView:numOfFFT/averageIncoherent:length(maxValue));

[value,pos]=max(aux(2,:));
valueCenteredForSavingPeak=aux(:,pos);
%%
fclose(fileID);
clear all;
