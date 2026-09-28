clear;
clc;
close all;

c = 299792458; %m/s , speed of light

%% Question 1 PRN 14

%Read in rinex file
fileinfo = rinexinfo('NIST00USA_R_20262310000_01D_30S_MO.rnx');
GPS_Data = rinexread('NIST00USA_R_20262310000_01D_30S_MO.rnx');

dataPRN14 = GPS_Data.GPS(GPS_Data.GPS.SatelliteID == 14, :); %Filter for PRN 5 data
timestep = seconds(30);
%Filtering, fill in spots that satellites isnt visible with NANS to make
%plots more accurate
dataPRN5 = retime(dataPRN5,'regular', 'fillwithmissing', 'TimeStep', timestep);

