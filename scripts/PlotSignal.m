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


% -----------------------
% FOR PLOTING SIGNALS
% -----------------------
%%
%OPEN FILE
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\s_6_d\PeaksIteration750.bin','r');
% fileID = fopen('prn_L1CA_32.bin','r');
%%
%READING COMPLEX SIGNALS
A_r = fread(fileID,[2 1000000000000000],'float32');  % Per senyals complexes
A_r = A_r(1,:)+1j*A_r(2,:);

%%
%READING REAL SIGNALS
% A = fread(fileID,10000000,'float32');  % Per senyals reals
% % A=A/max(A);

%%
%PLOTING COMPLEX SIGNALS
 figure
plot(1:length(A_r),abs(A_r));

%%
%PLOTING REAL
% hold on
% figure
% plot(A);
% % axis([0 10000 -1.5 1.5]);

%%
%CLOSE FILE
fclose(fileID);
% clear all;
%%
title ("50 Absolut peaks of real data");
% legend ("Average of 10 coherents");
%%
B = reshape(A_r,[311 1000]);
figure, plot(abs(B).^2)
% title ("50 Choerent Peaks, absolut value");
%%
C = sum(abs(B').^2);
C=C/max(C)
figure, plot(C)
title ("Sum of 50 Real Coherent Peaks");
%%
clear all;
clc;
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\s_6_d\PeaksIteration750.bin','r');
A_r = fread(fileID,[2 1000000000000000],'float32');  % Per senyals complexes
fclose(fileID);
A_r = A_r(1,:)+1j*A_r(2,:);
B = reshape(A_r,[311 1000]);
% C = sum(abs(B').^2);
% C=C/max(C);
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\s_6_d\Maximums750.bin','r');
Amax = fread(fileID,[5 1000000000000000],'single');
fclose(fileID);
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\s_6_r\Maximums750.bin','r');
Amax_r = fread(fileID,[5 1000000000000000],'single');
fclose(fileID);

% figure, plot(abs(B(:,300)))
% hold on;
fileID = fopen('C:\GNSS-REFLECTOMETRY-PROCESSING\results\s_6_r\PeaksIteration750.bin','r');
A_r = fread(fileID,[2 1000000000000000],'float32');  % Per senyals complexes
fclose(fileID);
A_r = A_r(1,:)+1j*A_r(2,:);
B_r = reshape(A_r,[311 1000]);
%%
mostres=floor(311/2);
temps=Amax(1)-mostres:1:Amax(1)+mostres;
temps_r=Amax_r(1)-mostres:1:Amax_r(1)+mostres;
temps=temps/(32*1.023e6)*(1023/1e-3);
temps_r=temps_r/(32*1.023e6)*(1023/1e-3);


% for i=1:250
mili=300;
directe=abs(B(:,mili));
reflectida=abs(B_r(:,mili));
factor_ = max(reflectida.^2);

figure,
plot(temps,(directe.^2)/factor_)
hold on;
plot(temps_r,(reflectida.^2)/factor_)


% C = sum(abs(B').^2);
% C_r = sum(abs(B_r').^2);
% 
% plot(temps,C/max(C))
% hold on;
% plot(temps_r,C_r/max(C_r))
% uiwait(msgbox('algunamerda'))
% end
grid minor
set(gca,'Fontsize',20)
title('Coherent power waveforms - GPS L1 C/A','Fontsize',20)
legend('Direct','Reflected','Location','NorthWest','Fontsize',20)
xlabel('Delay [C/A chips]','Fontsize',20)
ylabel('Normalized uncal. power','Fontsize',20)


%%
C = sum(abs(B').^2);
C=C/max(C);
plot(C)
