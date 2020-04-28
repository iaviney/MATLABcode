%2020 HW#1 Matlab and Plate Motions
% I Viney: 14 Feb 2020

clear; %clear all variables

%Create output file
outputfile=fopen('HW1_VINEY.txt', 'wt');

%Print header for output file
fprintf(outputfile,'HW1: Geos 419/519 Spring 2020: I. Viney\n'); 
fprintf(outputfile,'14 Feb 2020 \n');

% Read in Turcotte & Schubert (T&S) Plate Pole Positions & Rates (Table 1.6)
load Plate_Motions.mat ;
latd  = Velocities(1:30,1) ;
lond  = Velocities(1:30,2) ;
rate  = Velocities(1:30,3) ;

% Conv will convert degrees to radians for the trig functions
conv = (4*atan(1))/180. ;
lat = latd * conv ;
lon = lond * conv ;

% Convert all Poles in (lat, lon, rate) to Cartesian vectors (wx, wy, wz)
for i=1:30
wx(i)=rate(i)*cos(lon(i))*cos(lat(i));
wy(i)=rate(i)*sin(lon(i))*cos(lat(i));
wz(i)=rate(i)*sin(lat(i));
end

%Example 1
% Find the EU-NA relative motion pole and compare to the table
% EU-NA = -(AF-EU) + (AF-NA)

% AF-EU is plate pair #3 in the table, AF-NA is plate pair #2 in the table
EU_NA_x = -wx(3) + wx(2) ; % wx component of EU_NA pole
EU_NA_y = -wy(3) + wy(2) ; % wy component of EU_NA pole
EU_NA_z = -wz(3) + wz(2) ; % wz component of EU_NA pole

% Calculate rate as length of (wx,wy,wz) vector
rateEU_NA = sqrt(EU_NA_x^2 + EU_NA_y^2 + EU_NA_z^2) ; % rate, in deg/m.y.

% Convert back to Lat, Lon components from the cartesian using "Useful Equations.pdf"
latEU_NA = asin(EU_NA_z/rateEU_NA)/conv ; % latitude of pole (deg) (See "Useful Eqns")
lonEU_NA = atan2(EU_NA_y,EU_NA_x)/conv ;  % longitude of pole (deg)

% Print our EU-NA results to output file. Also print the EU-NA value from the table to compare:
fprintf(outputfile,'\n Example 1:\n');
fprintf(outputfile,'The EU-NA Pole is given by: \n');
fprintf(outputfile,'Lat(N) Lon(E)   Rate (deg/my):\n');
fprintf(outputfile,'%6.1f %6.1f %7.3f   based on EU-NA = -(AF-EU)+(AF-NA)\n',latEU_NA, lonEU_NA, rateEU_NA);
fprintf(outputfile,'%6.1f %6.1f %7.3f   from T&S Table 1.6 \n',latd(1), lond(1), rate(1)); % EU_NA is pole #1 in the table

%%%
% Question 3a
% Find the NA-PA relative motion pole and compare to the table
%%%

% Find the NA-PA relative motion pole and compare to the table
% NA-PA = -(EU-NA) + (EU-PA)

% EU-NA is plate pair #1 in the table, EU-PA is plate pair #18 in the table
NA_PA_x = -wx(1) + wx(18) ; % wx component of NA_PA pole
NA_PA_y = -wy(1) + wy(18) ; % wy component of NA_PA pole
NA_PA_z = -wz(1) + wz(18) ; % wz component of NA_PA pole

% Calculate rate as length of (wx,wy,wz) vector
rateNA_PA = sqrt(NA_PA_x^2 + NA_PA_y^2 + NA_PA_z^2) ; % rate, in deg/m.y.

% Convert back to Lat, Lon components from the cartesian using "Useful Equations.pdf"
latNA_PA = asin(NA_PA_z/rateNA_PA)/conv ; % latitude of pole (deg) (See "Useful Eqns")
lonNA_PA = atan2(NA_PA_y,NA_PA_x)/conv ;  % longitude of pole (deg)

% Print our NA-PA results to output file. Also print the NA-PA value from the table to compare:
fprintf(outputfile,'\n Question 3a:\n');
fprintf(outputfile,'The NA-PA Pole is given by: \n');
fprintf(outputfile,'Lat(N) Lon(E)   Rate (deg/my):\n');
fprintf(outputfile,'%6.1f %6.1f %7.3f   based on NA-PA = -(EU-NA)+(EU-PA)\n',latNA_PA, lonNA_PA, rateNA_PA);
fprintf(outputfile,'%6.1f %6.1f %7.3f   from T&S Table 1.6 \n',latd(9), lond(9), rate(9)); % NA_PA is pole #9 in the table

%%%
% Question 3b
% Find the PA-PA relative motion pole using a circuit with at least 2 other plates
%%%

% PA-PA = (NA-PA) + (AF-NA) + (AU-AF) + (PA-AU)

% NA-PA is plate pair #9 in the table, AF-NA is plate pair #2 in the table,
% AU-AF is plate pair #23 in the table, PA-AU is plate pair #17 in the
% table
PA_PA_x = wx(9) + wx(2) + wx(23) + wx(17) ; % wx component of PA_PA pole
PA_PA_y = wy(9) + wy(2) + wy(23) + wy(17) ; % wy component of PA_PA pole
PA_PA_z = wz(9) + wz(2) + wz(23) + wz(17) ; % wz component of PA_PA pole

% Calculate rate as length of (wx,wy,wz) vector
ratePA_PA = sqrt(PA_PA_x^2 + PA_PA_y^2 + PA_PA_z^2) ; % rate, in deg/m.y.

% Convert back to Lat, Lon components from the cartesian using "Useful Equations.pdf"
latPA_PA = asin(PA_PA_z/ratePA_PA)/conv ; % latitude of pole (deg) (See "Useful Eqns")
lonPA_PA = atan2(PA_PA_y,PA_PA_x)/conv ;  % longitude of pole (deg)

% Print our PA-PA results to output file
fprintf(outputfile,'\n Question 3b:\n');
fprintf(outputfile,'The PA-PA Pole is given by: \n');
fprintf(outputfile,'Lat(N) Lon(E)   Rate (deg/my):\n');
fprintf(outputfile,'%6.1f %6.1f %7.3f   based on PA-PA = (NA-PA) + (AF-NA) + (AU-AF) + (PA-AU)\n',latPA_PA, lonPA_PA, ratePA_PA);

%%%
% Question 3c
% Find the AF-PA relative motion pole 
%%%

% AF-PA = (AF-EU) + (EU-PA)

% AF-EU is plate pair #3 in the table, EU-PA is plate pair #18 in the table
AF_PA_x = wx(3) + wx(18) ; % wx component of AF_PA pole
AF_PA_y = wy(3) + wy(18) ; % wy component of AF_PA pole
AF_PA_z = wz(3) + wz(18) ; % wz component of AF_PA pole

% Calculate rate as length of (wx,wy,wz) vector
rateAF_PA = sqrt(AF_PA_x^2 + AF_PA_y^2 + AF_PA_z^2) ; % rate, in deg/m.y.

% Convert back to Lat, Lon components from the cartesian using "Useful Equations.pdf"
latAF_PA = asin(AF_PA_z/rateAF_PA)/conv ; % latitude of pole (deg) (See "Useful Eqns")
lonAF_PA = atan2(AF_PA_y,AF_PA_x)/conv ;  % longitude of pole (deg)

% Print our AF-PA results to output file
fprintf(outputfile,'\n Question 3c:\n');
fprintf(outputfile,'The AF-PA Pole is given by: \n');
fprintf(outputfile,'Lat(N) Lon(E)   Rate (deg/my):\n');
fprintf(outputfile,'%6.1f %6.1f %7.3f   based on AF-PA = (AF-EU) + (EU-PA)\n',latAF_PA, lonAF_PA, rateAF_PA);

fclose(outputfile); % This closes file outputfile, a necessary last step. Use only once!

