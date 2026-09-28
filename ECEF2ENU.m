function C_ECEF2ENU = ECEF2ENU(ref_lat_deg, ref_lon_deg)
%ECEF2ENU builds the DCM matrix to rotate vectors from ECEF frame to ENU
%frame, it just uses the known DCM matrix provided in lecture

%Assigning variables to make matrix look less cluttered
slamb= sin(deg2rad(ref_lon_deg));
clamb = cos(deg2rad(ref_lon_deg));

sphi = sin(deg2rad(ref_lat_deg));
cphi = cos(deg2rad(ref_lat_deg));

C_ECEF2ENU = [-slamb clamb 0; -sphi*clamb -sphi*slamb cphi; cphi*clamb cphi*slamb sphi];

end