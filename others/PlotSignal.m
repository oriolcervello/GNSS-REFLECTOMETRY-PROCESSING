% -----------------------
% FOR PLOTING SIGNALS
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
%READING REAL SIGNALS
A = fread(fileID,10000000,'float32');  % Per senyals reals
% A=A/max(A);

%%
%PLOTING COMPLEX SIGNALS
 hold on
% figure
plot(1:length(A),abs(A));

%%
%PLOTING REAL
% hold on
figure
plot(A);
% axis([0 10000 -1.5 1.5]);

%%
%CLOSE FILE
fclose(fileID);
clear all;
%%
title ("50 Absolut peaks of real data");
% legend ("Average of 10 coherents");
%%
 B = reshape(A,[311 80]);
figure, plot(abs(B').^2)
title ("50 Choerent Peaks, absolut value");
%%
C = sum(abs(B').^2);
figure, plot(C)
title ("Sum of 50 Real Coherent Peaks");



