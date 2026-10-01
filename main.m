clear;
clc;
close all;

c = 299792458; %m/s , speed of light
NIST_ECEF = [-1288398.567 -4721696.932	4078625.350]; %Nist ecef location

%% Question 1 PRN 14

%Read in rinex file
fileinfo = rinexinfo('NIST00USA_R_20262310000_01D_30S_MO.rnx');
GPS_Data = rinexread('NIST00USA_R_20262310000_01D_30S_MO.rnx');
%Read ephemeris values
[gps_ephem,ionoparams] = read_clean_GPSbroadcast('brdc2310.26n');

dataPRN14 = GPS_Data.GPS(GPS_Data.GPS.SatelliteID == 14, :); %Filter for PRN 14 data
timestep = seconds(30);
%Filtering, fill in spots that satellites isnt visible with NANS to make
%plots more accurate
dataPRN14 = retime(dataPRN14,'regular', 'fillwithmissing', 'TimeStep', timestep);

%Need to manipulate RINEX data so we can use findExpectedRange function
% Define the gps start time
gpsStart14 = datetime(1980, 1, 6, 0, 0, 0);

% Extract timestamps 
timestamps14 = dataPRN14.Time; 

% Compute elapsed days and divide by 7
daysSinceStart14 = days(timestamps14 - gpsStart14);
gps_week14 = floor(daysSinceStart14 / 7);

% Determine day of week
dayNum14 = day(dataPRN14.Time, 'dayofweek') - 1; % -1 is so sunday =0

% Convert to TOW
receiverTime14 = dayNum14 * 86400 + seconds(timeofday(dataPRN14.Time));

expectedRangefromRinex14 = findExpectedRange(NIST_ECEF,receiverTime14,gps_week14,gps_ephem,14);

expectedRangefromRinex14(isnan(dataPRN14.C1C)) = NaN;

dPR0 = dataPRN14.C1C - expectedRangefromRinex14; %find error between c1c and expected range
fprintf('First dPR0 value: %.4f m\n', dPR0(1));
fprintf('Last dPR0 value:  %.4f m\n', dPR0(end));


figure('Name','prn26pseudoandexpected');

plot(receiverTime14./3600 - 72, dataPRN14.C1C);
hold on;
plot(receiverTime14./3600 - 72, expectedRangefromRinex14,'r--')
title('PRN 14 C1C pseudorange and Expected Range')
xlabel('Time (hr)')
ylabel('Range (m)')
legend('C1C Pseudorange', 'Expected Range')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure('Name','prn14pseudoexpecteddiff');
plot(receiverTime14./3600 - 72,dPR0)
title('Error between PRN 14 C1C and Expected Range')
xlabel('Time (hr)')
ylabel('Error (m)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%% Question 2

[health14,satPos14,satVel14,satClkCorr14,relCorr14,junk14,tgd14] = eph2pvt2025(gps_ephem,[gps_week14 receiverTime14],14);

figure("Name","Satellite 14 Clock Correction Vs Time")
plot(receiverTime14./3600 - 72,satClkCorr14);
title('Satellite 14 Clock Bias Vs Time')
xlabel('Time (hr)')
ylabel('Satellite Clock Bias (m)')

dPR1 = dataPRN14.C1C - (expectedRangefromRinex14 - satClkCorr14);
fprintf('First dPR1 value: %.4f m\n', dPR1(1));
fprintf('Last dPR1 value:  %.4f m\n', dPR1(end));

figure("Name","Satellite 14 with simple Clock Correction")
plot(receiverTime14./3600 - 72,dPR1);
title('Error between PRN 14 and Expected Range with simple clock correction applied')
xlabel('Time (hr)')
ylabel('Error (m)')

%% Question 3 Updated eph2pvt2025 function to include relativistic correction

dPR2 = dataPRN14.C1C - (expectedRangefromRinex14 - satClkCorr14 - relCorr14);
fprintf('First dPR2 value: %.4f m\n', dPR2(1));
fprintf('Last dPR2 value:  %.4f m\n', dPR2(end));

figure("Name","Relativistic Clock Correction")
plot(receiverTime14./3600 - 72,relCorr14);
title('Relativistic Clock Correction PRN 14')
xlabel('Time (hr)')
ylabel('Relativistic Clock Correction (m)')

figure("Name","Satellite 14 with simple, relativistic Clock Correction")
plot(receiverTime14./3600 - 72,dPR2);
title('Error between PRN 14 and Expected Range with simple and relativistic clock correction applied')
xlabel('Time (hr)')
ylabel('Error (m)')

%% Question 4

%find elevation of satellite
[AZ14, EL14, RANGE14] = compute_azelrange(NIST_ECEF, satPos14);

zd_NIST = 2;

tropo = tropomodel(zd_NIST,EL14);

dPR3 = dataPRN14.C1C - (expectedRangefromRinex14 - satClkCorr14 - relCorr14 + tropo);
fprintf('First dPR3 value: %.4f m\n', dPR3(1));
fprintf('Last dPR3 value:  %.4f m\n', dPR3(end));

figure("Name","troposphere Clock Correction")
plot(receiverTime14./3600 - 72,tropo);
title('Tropo clock correction')
xlabel('Time (hr)')
ylabel('Tropo correction (m)')


figure("Name","Satellite 14 with simple, relativistic, and troposphere Clock Correction")
plot(receiverTime14./3600 - 72,dPR3);
title('Error between PRN 14 and Expected Range with simple, relativistic, tropo clock correction applied')
xlabel('Time (hr)')
ylabel('Error (m)')


%Option problem 3 work

NIST_LLA = ecef2lla(NIST_ECEF);
NIST_H = 1658.8; %Orthometric height of nist

for i = 1:length(EL14)

[tropoUNB(i,1), ~ ,~ ,~ ,~]= UNB3M(rad2deg(NIST_LLA(1)),NIST_H,231,rad2deg(EL14(i)));

end

dPR3UNB = dataPRN14.C1C - (expectedRangefromRinex14 - satClkCorr14 - relCorr14 + tropoUNB);


figure("Name","Troposphere Clock Correction simple vs UNB3")
plot(receiverTime14./3600 - 72,tropo); hold on;
plot(receiverTime14./3600 - 72, tropoUNB)
title('Tropo clock correction')
xlabel('Time (hr)')
ylabel('Tropo correction (m)')
legend('Simple Model Eq 5.43', 'UNB Model')


figure("Name","Troposphere UNB3 correction applied")
plot(receiverTime14./3600 - 72, dPR3UNB)
title('Tropo clock correction')
xlabel('Time (hr)')
ylabel('Tropo correction (m)')

%% Question 5


[PRIF12, iono14] = ionocorr(dataPRN14.C1C, 1575.42e6, dataPRN14.C2L, 1227.60e6);

dPR4 = PRIF12 - (expectedRangefromRinex14 - satClkCorr14 - relCorr14 + tropo);
fprintf('First dPR4 value: %.4f m\n', dPR4(1));
fprintf('Last dPR4 value:  %.4f m\n', dPR4(end));


figure("Name","Ionospheric correction (m)")
plot(receiverTime14./3600 - 72, iono14)
title('Iono clock correction')
xlabel('Time (hr)')
ylabel('Iono correction (m)')


figure("Name","Ionospheric correction applied")
plot(receiverTime14./3600 - 72, dPR4)
title('Iono clock correction applied (dPR4)')
xlabel('Time (hr)')
ylabel('Satellite Range Difference (m)')


%% Question 6 
figure("Name","Comparing Different Correction")
plot(receiverTime14./3600 - 72, dPR1)
hold on
plot(receiverTime14./3600 - 72, dPR2)
plot(receiverTime14./3600 - 72, dPR3)
%plot(receiverTime14./3600 - 72, dPR3UNB)
plot(receiverTime14./3600 - 72, dPR4)
title('Iono clock correction applied (dPR4')
xlabel('Time (hr)')
ylabel('Satellite Range Difference (m)')
legend('dpr1','dpr2','dpr3','dpr4')

%% Question 7

 [MP1,CMC1] = mpath(dataPRN14.C1C , dataPRN14.L1C, 1575.42e6, dataPRN14.L2W, 1227.60e6); %C1C call
 [MP2,CMC2] = mpath(dataPRN14.C2W , dataPRN14.L2W, 1227.60e6, dataPRN14.L1C, 1575.42e6); %C1W call, need a different second frequency
 [MP3,CMC3] = mpath(dataPRN14.C2L , dataPRN14.L2L, 1227.60e6, dataPRN14.L1C, 1575.42e6); %C2L call
 [MP4,CMC4] = mpath(dataPRN14.C5Q , dataPRN14.L5Q, 1176.45e6, dataPRN14.L2W, 1227.60e6); %C5Q call



figure("Name","CMC of C1C")

subplot(2,1,1)
plot(receiverTime14./3600 - 72, CMC1)
title('CMC1')
xlabel('Time (hr)')
ylabel('Satellite Range Difference (m)')

subplot(2,1,2)
plot(receiverTime14./3600 - 72, dataPRN14.S1C)
title('S1C')
xlabel('Time (hr)')
ylabel('SNR (dbHz)')

figure("Name","CMC1 of C2W")

subplot(2,1,1)
plot(receiverTime14./3600 - 72, CMC2)
title('CMC of C2W')
xlabel('Time (hr)')
ylabel('Satellite Range Difference (m)')

subplot(2,1,2)
plot(receiverTime14./3600 - 72, dataPRN14.S2W)
title('S2W')
xlabel('Time (hr)')
ylabel('SNR (dbHz)')

figure("Name","CMC of C2L")
subplot(2,1,1)
plot(receiverTime14./3600 - 72, CMC3)
title('CMC of C2L')
xlabel('Time (hr)')
ylabel('Satellite Range Difference (m)')

subplot(2,1,2)
plot(receiverTime14./3600 - 72, dataPRN14.S1L)
title('S2L')
xlabel('Time (hr)')
ylabel('SNR (dbHz)')

figure("Name","CMC of C5Q")

subplot(2,1,1)
plot(receiverTime14./3600 - 72, CMC4)
title('CMC of C5Q')
xlabel('Time (hr)')
ylabel('Satellite Range Difference (m)')

subplot(2,1,2)
plot(receiverTime14./3600 - 72, dataPRN14.S5Q)
title('S5Q')
xlabel('Time (hr)')
ylabel('SNR (dbHz)')