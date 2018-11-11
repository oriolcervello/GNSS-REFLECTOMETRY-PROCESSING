fileID = fopen('C:\Users\Documents\GNSS-REFLECTOMETRY-PROCESSING\results\PeaksIteration2.bin','r');
% fileID = fopen('PeaksIteration2.bin','r');
%%
%Reading comlex signal
A = fread(fileID,[2 10000000],'float32');  % Per senyals complexes

A = A(1,:)+1j*A(2,:);


%%
%Complex signals
numofFFTs=20;
peakrangetosave=311;
% hold on
figure
plot((1:peakrangetosave*numofFFTs),A(1,1:peakrangetosave*numofFFTs));
fclose(fileID);

%%
fclose(fileID);

