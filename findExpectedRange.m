function [expectedRange] = findExpectedRange(RecieverECEF, Tr, gpsWeek, gps_ephem,prn);
%{
Author: Yoan Serafimov

INPUTS:
RecieverECEF - ECEF positions of receiver
Tr = reciever time in TOW
gps_ephem = ephemris matrix for all gps satellites
prn - prn number for the satellite you want the expected range for
gpsWeek - gps week number

OUTPUTS:
expected range from transmission time through algorithm explained in hw3


%}



%Filter gps_ephem matrix for chosen PRN
gps_ephem = gps_ephem(gps_ephem(:,1) == prn, :);



w_e = 7.2921150e-5; % earths rotation rate rad/s\
c = 299792458; %m/s




%Compute PRN  position with ephemeris at times in the sp3 file
[~,satAt_Tr,~,~,~,~] = eph2pvt2025(gps_ephem,[gpsWeek Tr],prn);

%Pre-allocation for speed (very computationally expensive problem)!
N = length(Tr);
expectedRange = zeros(N,1);

Rmatrix = @(phi) [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];


for i=1:length(Tr) %compute expected range iteratively
    expectedRange(i) = norm(satAt_Tr(i,:) - RecieverECEF);
    
    for j = 1:3 %repeat algorithm 3 times
    %Step 3
    Tt = Tr(i) - expectedRange(i)./c;

    %Step 4
    [~,satAt_Tt,~,~,~,~] = eph2pvt2025(gps_ephem,[gpsWeek(i) Tt],prn);

    %step 5

    phi = w_e.*(Tr(i) - Tt);

    satAt_Tr(i,:) = (Rmatrix(phi)*(satAt_Tt)')';

    %step 6
    expectedRange(i) = norm(satAt_Tr(i,:) - RecieverECEF);
    end
    
end


end

