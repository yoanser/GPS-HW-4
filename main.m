clear;
clc;
close all;

c = 299792458; %m/s , speed of light

%% Question 1 PRN 05

%Read in rinex file
fileinfo = rinexinfo('NIST00USA_R_20262310000_01D_30S_MO.rnx');
GPS_Data = rinexread('NIST00USA_R_20262310000_01D_30S_MO.rnx');

dataPRN5 = GPS_Data.GPS(GPS_Data.GPS.SatelliteID == 5, :); %Filter for PRN 5 data
timestep = seconds(30);
%Filtering, fill in spots that satellites isnt visible with NANS to make
%plots more accurate
dataPRN5 = retime(dataPRN5,'regular', 'fillwithmissing', 'TimeStep', timestep);

%Plot! Chose L1C for carrier phase,
figure("Name","PRN5 Pseudorange");
plot(dataPRN5,"Time",'C1C')
ylabel('C1C Pseudorange (m)')


cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure("Name","PRN5 Signal to Noise Ratio")
plot(dataPRN5,"Time","S1C")
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure("Name","PRN5 L1C")
L1C_Lambda = c/(1575.42e6); %Find wavelength of L1 band
dataPRN5_L1CinMeters = dataPRN5;
dataPRN5_L1CinMeters.L1C = dataPRN5_L1CinMeters.L1C .* L1C_Lambda; % Converting L1C carrier phase to meters
plot(dataPRN5_L1CinMeters,"Time","L1C")
ylabel("L1C (meters) ")
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);


figure("Name","PRN5 C1W")
plot(dataPRN5,"Time","C1W");
ylabel("C1W (meters)");
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%% Question 2

%Read ephemeris values
[gps_ephem,ionoparams] = read_clean_GPSbroadcast('brdc2310.26n');

%Read data
sp3=read_sp3('IGS0OPSFIN_20262310000_01D_15M_ORB.SP3');

%parse date for prn5 only
prn5sp3data = sp3(sp3(:,3)==5 , :);

time_in_hours = prn5sp3data(:,2)./(3600) - 72;%Make a time array in hours from start of day;

%Compute PRN 5 position with ephemeris at times in the sp3 file
[health5,satPos5,satVel5,satClkCorr5,junk5,tgd5] = eph2pvt2025(gps_ephem,prn5sp3data(:,1:2),5);

%Plot exact PRN 5 position

figure('Name','Exact PRN 5 Position overlaid with ephemeris')

subplot(3,1,1);
plot(time_in_hours,prn5sp3data(:,4).*1000);
title('GPS 5 exact position vs time of day')
hold on;
plot(time_in_hours,satPos5(:,1),'r--');
xlabel('time (hr)')
ylabel('X position (m)')
legend('Exact position', 'Ephemeris-based position')

subplot(3,1,2)
plot(time_in_hours,prn5sp3data(:,5).*1000,'g');
hold on;
plot(time_in_hours,satPos5(:,2),'r--');
xlabel('time (hr)')
ylabel('Y position (m)')
legend('Exact position', 'Ephemeris-based position')


subplot(3,1,3)
plot(time_in_hours,prn5sp3data(:,6).*1000,'k');
hold on;
plot(time_in_hours,satPos5(:,3),'r--');
xlabel('time (hr)')
ylabel('Z position (m)')
legend('Exact position', 'Ephemeris-based position')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%Plot ephemeris based satellite position
figure('Name','Ephemeris Based PRN 5 Position')

subplot(3,1,1);
plot(time_in_hours,satPos5(:,1));
title('GPS 5 Ephemeris Based Position')
xlabel('time (hr)')
ylabel('X position (m)')

subplot(3,1,2)
plot(time_in_hours,satPos5(:,2),'r');
xlabel('time (hr)')
ylabel('Y position (m)')


subplot(3,1,3)
plot(time_in_hours,satPos5(:,3),'k');
xlabel('time (hr)')
ylabel('Z position (m)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure('Name','Difference Between Ephem and Exact')

subplot(3,1,1);
plot(time_in_hours,prn5sp3data(:,4).*1000 - satPos5(:,1));
title('GPS 5 position error between ephemeris and exact vs time of day')
xlabel('time (hr)')
ylabel('X position (m)')

subplot(3,1,2)
plot(time_in_hours,prn5sp3data(:,5).*1000 - satPos5(:,2),'r');
xlabel('time (hr)')
ylabel('Y position (m)')


subplot(3,1,3)
plot(time_in_hours,prn5sp3data(:,6).*1000 - satPos5(:,3),'k');
xlabel('time (hr)')
ylabel('Z position (m)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%part c find the clock bias of PRN 5
%Filter gps_ephem matrix for PRN 5
gps_ephem5 = gps_ephem(gps_ephem(:,1) == 5, :);

%Extract toc and toe values
tobs = prn5sp3data(:,2);
toc = gps_ephem5(:,20);
a0 = gps_ephem5(:,21);
a1 = gps_ephem5(:,22);

time_in_hoursephem = gps_ephem5(:,17)./(3600);%Make a time array in hours from start of day;

%calculate satelitte clock bias
for i = 1:length(tobs) %sample throuhg sp3 time vector
currentTime = tobs(i);
[~, minIdx] = min(abs(toc - currentTime)); %Choose indexes that correspond most closely to current time

clockBias5(i) = a0(minIdx) + a1(minIdx).*(currentTime - toc(minIdx)); %find clock bias

end



figure('Name','Clock Bias in Time PRN 5')
plot(time_in_hours,clockBias5)
title('PRN 5 Clock Bias')
xlabel("Time (hours)")
ylabel("Satellite Clock Bias (Seconds)")
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%Find clock bias with relativistic correction
%find dot product between r and v:

dotRVprn5 = dot(satPos5,satVel5,2) %These position and velocity vectors are from the ephemeris broadcast

%Find relativistic correction!

c = 299792458; %m/s
relativisticCorrection = (2.*dotRVprn5)./c^2;

figure('Name','Relativistic Time Correction')
plot(time_in_hours,relativisticCorrection)
title('Relativist Time Correction for PRN 5')
ylabel('Seconds')
xlabel('Time (hr)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%% Question 3
c = 299792458; %m/s
NIST_ECEF = [-1288398.567 -4721696.932	4078625.350]; %Nist ecef location

[AZ5, EL5, RANGE5] = compute_azelrange(NIST_ECEF, satPos5); %satPos5 was found in problem 2

figure('Name','prn5rangeandelevation')
subplot(2,1,1)
plot(time_in_hours,EL5,'b'); 
title('Plotting elevation and Range of PRN5 based on ephemeris')
ylabel('El (deg)')
xlabel('time (hr)')


subplot(2,1,2)
plot(time_in_hours,RANGE5,'k'); 
ylabel('Range (m)')
xlabel('time (hr)') 
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

%% Find expected range and plot difference between ephemeris range and rotation range
expectedRange = findExpectedRange(NIST_ECEF,prn5sp3data(:,2),prn5sp3data(:,1),gps_ephem,5);



figure('Name','prn5rangeandelevationwexpected')
subplot(2,1,1)
plot(time_in_hours,EL5,'b'); 
title('Plotting elevation and Range of PRN5 based on ephemeris')
ylabel('El (deg)')
xlabel('time (hr)')


subplot(2,1,2)
plot(time_in_hours,RANGE5,'k'); 
hold on;
plot(time_in_hours,expectedRange,'r--')
ylabel('Range (m)')
xlabel('time (hr)') 
legend('Range from sp3', 'Expected Range')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure('Name','errorEphemandExpected')
plot(time_in_hours,(RANGE5' - expectedRange))
ylabel('Error (m)')
xlabel('time (hr)')
title('Error between Ephemeris Range and expected range')
maxDiff5 = max(abs(RANGE5' - expectedRange))
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);


%% Question 4
%Need to manipulate RINEX data so we can use findExpectedRange function
% Define the gps start time
gpsStart5 = datetime(1980, 1, 6, 0, 0, 0);

% Extract timestamps 
timestamps5 = dataPRN5.Time; 

% Compute elapsed days and divide by 7
daysSinceStart5 = days(timestamps5 - gpsStart5);
gps_week5 = floor(daysSinceStart5 / 7);

% Determine day of week
dayNum5 = day(dataPRN5.Time, 'dayofweek') - 1; % -1 is so sunday =0

% Convert to TOW
receiverTime5 = dayNum5 * 86400 + seconds(timeofday(dataPRN5.Time));

expectedRangefromRinex5 = findExpectedRange(NIST_ECEF,receiverTime5,gps_week5,gps_ephem,5);

expectedRangefromRinex5(isnan(dataPRN5.C1C)) = NaN;

figure('Name','prn5pseudoandexpected');

plot(receiverTime5./3600 - 72, dataPRN5.C1C);
title('PRN 5 C1C pseudorange and Expected Range')
hold on;
plot(receiverTime5./3600 - 72, expectedRangefromRinex5,'r--')

xlabel('Time (hr)')
ylabel('Range (m)')
legend('C1C Pseudorange', 'Expected Range')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure('Name','prn5pseudoandexpecteddiff');
plot(receiverTime5./3600 - 72,dataPRN5.C1C - expectedRangefromRinex5)
title('Error between PRN 5 C1C and Expected Range')
xlabel('Time (hr)')
ylabel('Error (m)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);


% Test my theory with another satellite (PRN 13) to see if the reciever
% clock bias is the main source of error between pseudorange and expected
% range


dataPRN7 = GPS_Data.GPS(GPS_Data.GPS.SatelliteID == 7, :); %Filter for PRN 7 data

timestep = seconds(30);
%Filtering, fill in spots that satellites isnt visible with NANS to make
%plots more accurate
dataPRN7 = retime(dataPRN7,'regular', 'fillwithmissing', 'TimeStep', timestep);

%Need to manipulate RINEX data so we can use findExpectedRange function
% Define the gps start time
gpsStart7 = datetime(1980, 1, 6, 0, 0, 0);

% Extract timestamps 
timestamps7 = dataPRN7.Time; 
% Compute elapsed days and divide by 7
daysSinceStart7 = days(timestamps7 - gpsStart7);
gps_week7 = floor(daysSinceStart7 / 7);

% Determine day of week
dayNum7 = day(dataPRN7.Time, 'dayofweek') - 1; % -1 is so sunday =0

% Convert to TOW
receiverTime7 = dayNum7 * 86400 + seconds(timeofday(dataPRN7.Time));

expectedRangefromRinex7 = findExpectedRange(NIST_ECEF,receiverTime7,gps_week7,gps_ephem,7);

expectedRangefromRinex7(isnan(dataPRN7.C1C)) = NaN;

figure('Name','prn7pseudoandexpected');

plot(receiverTime7./3600 - 72, dataPRN7.C1C);
hold on;
plot(receiverTime7./3600 - 72, expectedRangefromRinex7,'r--')
title('PRN 7 C1C pseudorange and Expected Range')
xlabel('Time (hr)')
ylabel('Range (m)')
legend('C1C Pseudorange', 'Expected Range')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure('Name','prn7pseudoandexpecteddiff');
plot(receiverTime7./3600 - 72,dataPRN7.C1C - expectedRangefromRinex7)
title('Error between PRN 7 C1C and Expected Range')
xlabel('Time (hr)')
ylabel('Error (m)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);




%% Question 5 let's repeat question 4 for PRN 26


dataPRN26 = GPS_Data.GPS(GPS_Data.GPS.SatelliteID == 26, :); %Filter for PRN 26 data

timestep = seconds(30);
%Filtering, fill in spots that satellites isnt visible with NANS to make
%plots more accurate
dataPRN26 = retime(dataPRN26,'regular', 'fillwithmissing', 'TimeStep', timestep);

%Need to manipulate RINEX data so we can use findExpectedRange function
% Define the gps start time
gpsStart7 = datetime(1980, 1, 6, 0, 0, 0);

% Extract timestamps 
timestamps26 = dataPRN26.Time; 

% Compute elapsed days and divide by 7
daysSinceStart26 = days(timestamps26 - gpsStart7);
gps_week26 = floor(daysSinceStart26 / 7);

% Determine day of week
dayNum26 = day(dataPRN26.Time, 'dayofweek') - 1; % -1 is so sunday =0

% Convert to TOW
receiverTime26 = dayNum26 * 86400 + seconds(timeofday(dataPRN26.Time));

expectedRangefromRinex26 = findExpectedRange(NIST_ECEF,receiverTime26,gps_week26,gps_ephem,26);

expectedRangefromRinex26(isnan(dataPRN26.C1C)) = NaN;

figure('Name','prn26pseudoandexpected');

plot(receiverTime26./3600 - 72, dataPRN26.C1C);
hold on;
plot(receiverTime26./3600 - 72, expectedRangefromRinex26,'r--')
title('PRN 26 C1C pseudorange and Expected Range')
xlabel('Time (hr)')
ylabel('Range (m)')
legend('C1C Pseudorange', 'Expected Range')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);

figure('Name','prn26pseudoexpecteddiff');
plot(receiverTime26./3600 - 72,dataPRN26.C1C - expectedRangefromRinex26)
title('Error between PRN 26 C1C and Expected Range')
xlabel('Time (hr)')
ylabel('Error (m)')
cleanFileName = strrep(gcf().Name, ' ', '_');  
fileName = sprintf('%s.png', cleanFileName);
exportgraphics(gcf, fileName, 'Resolution', 300);