function [tropo] = tropomodel(zd,elevation)


tropo = zd.*(1./sqrt(1-(cosd(elevation)./1.001).^2)); %Implementing Eq ME 5.42, no flat surface approximation

tropo = tropo'; %makes tropo a column vector so its consistent with the other outputs

end