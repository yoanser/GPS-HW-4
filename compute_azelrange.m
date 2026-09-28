function [AZ, EL, RANGE] = compute_azelrange(userECEF, satECEF)

%For loop so we can handle array of satellite positions
for i=1:height(satECEF)

%Find range (It is simply the magnitude of the difference between the user
%and sat vectors)
vectorObstoSV = satECEF(i,:) - userECEF;
RANGE(i) = norm(vectorObstoSV);

%Find relative vector in ENU frame
ObsGeodetic = ecef2lla(userECEF) ;
C_ECEF2ENU = ECEF2ENU(ObsGeodetic(1), ObsGeodetic(2));

% Convert the vector to ENU coordinates
enuVector = C_ECEF2ENU * vectorObstoSV';

% Calculate azimuth and elevation angles (degrees)
AZ(i) = wrapTo360(atan2d(enuVector(1),enuVector(2))); %Use wrap to 360 to convert negative angels
EL(i) = asind(enuVector(3)./ norm(enuVector));

end

end