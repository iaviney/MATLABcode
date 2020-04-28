% 2020 HW#2 Matlab and Plate Tectonics 2
%
% I VINEY: 19 Feb 2020 

%%%
% Setup
%%% 

clear ; % clear all variables
close all ; % close any open files

% Open an output file. 
outputfile = fopen('HW2_VINEY.txt', 'wt');

% Print a header to output file
fprintf(outputfile,'HW2: Geos 419/519 Spring 2020: I. VINEY\n'); 
fprintf(outputfile,'20 Feb 2020 \n');

% Read in Turcotte & Schubert (T&S) Plate Pole Positions & Rates (Table 1.6)
load Plate_Motions.mat ; % This creates a 30x3 variable called Velocities with the T&S data
latTSd  = Velocities(1:30,1) ;  % latd is a 30x1 array of pole latitudes (in degrees)
lonTSd  = Velocities(1:30,2) ;  % lond is a 30x1 array of pole longitudes (in degrees)
rateTS  = Velocities(1:30,3) ;  % rate is a 30x1 array or plate rates (deg/m.y.)

% Conv will convert degrees to radians for the trig functions
conv = (4*atan(1))/180 ;

% Here we convert the table data in degrees to data in radians
latTSr = latTSd * conv ;  % lat is a 30x1 array of pole latitudes (in radians)
lonTSr = lonTSd * conv ;  % lon is a 30x1 array of pole longitudes (in radians)

% Convert all table data Poles in (lat, lon, rate) (in radians) to Cartesian vectors (wx, wy, wz)
for i = 1:30 % for the 30 plate-pair poles in T&S Table 1.6
    wxTS(i) = rateTS(i) * cos(lonTSr(i)) * cos(latTSr(i)) ; 
    wyTS(i) = rateTS(i) * sin(lonTSr(i)) * cos(latTSr(i)) ;
    wzTS(i) = rateTS(i) * sin(latTSr(i)) ; 
end

% Earth radius in km;
radius = 6371 ;  PI = 4*atan(1);
% scale converts deg to km
scale = (1./360.) * 2*(4*atan(1)) * radius ; 

% ABSOLUTE plate velocity pole for the Pacific plate in a fixed Hot Spot 
% Reference Frame (HS3-NUVEL1A):

% Latitude
latPA_HSd = -61.5453; %degrees
% Longitude
lonPA_HSd = 90.6174; %degrees
% Rate
ratePA_HS = 1.0633; %degree/my

% Haleakala plate location in spherical coordinates (in degrees)
latHd = 20.7 ;
lonHd = -156.3 ;
rateH = 1 ; % using  unit radius (rate information from pole rotation)


% Print PA_HS and H to HW2.out:
fprintf(outputfile,'\nThe absolute plate velocity of the PA_HS Pole: \n');
fprintf(outputfile,'Lat(N), Lon(E),   Rate (deg/my):\n');
fprintf(outputfile,'%1.1f  %5.1f  %9.3f\n',latPA_HSd, lonPA_HSd, ratePA_HS);
fprintf(outputfile,'\nThe location of Haleakala: \n');
fprintf(outputfile,'Lat(N), Lon(E),   Rate (unit rate):\n');
fprintf(outputfile,'%1.1f  %8.1f  %7.3f\n',latHd, lonHd, rateH);

% Setup complete

%%%
% Question 1a
% Calculate the cartesian velocity of Haleakala
%%%

% Converting the rotation pole spherical coordinates to radians:
latPA_HSr=latPA_HSd*conv;
lonPA_HSr=lonPA_HSd*conv;

% Converting the Haleakala point spherical cooridinates to radians:
latHr = latHd * conv ; 
lonHr = lonHd * conv ;

%Converting rotation pole coordinates from spherical (radians) into
%cartesian:
wxPA_HS = ratePA_HS * cos(lonPA_HSr) * cos(latPA_HSr) ; 
wyPA_HS = ratePA_HS * sin(lonPA_HSr) * cos(latPA_HSr) ;
wzPA_HS = ratePA_HS * sin(latPA_HSr) ; 

%Converting Haleakala point spherical cooridinates (radians) to cartesian: 
wxH = rateH * cos(lonHr) * cos(latHr) ; 
wyH = rateH * sin(lonHr) * cos(latHr) ;
wzH = rateH * sin(latHr) ; 

%Calculating absolute motion of Haleakala plate (degrees/my)
Hvx=(wyPA_HS*wzH)-(wzPA_HS*wyH);
Hvy=((-wxPA_HS)*wzH)+(wzPA_HS*wxH);
Hvz=(wxPA_HS*wyH)-(wyPA_HS*wxH);

%Scale from degrees/my to km/my (AKA mm/yr)
sHvx=scale*Hvx ;
sHvy=scale*Hvy ;
sHvz=scale*Hvz ;

%Calculating the length of the velocity 
lengthHv=sqrt(sHvx^2 + sHvy^2 + sHvz^2);

%Printing to answer file
fprintf(outputfile,'\nQuestion 1a:\n');
fprintf(outputfile,'\nHaleakala Cartesian Velocity:\n');
fprintf(outputfile,'vx,      vy,     vz   length:\n');
fprintf(outputfile,'%0.1f  %4.1f  %4.1f  %4.3f\n',sHvx, sHvy, sHvz, lengthHv);

%Testing my script using the North Pole:

% Converting the North Pole spherical cooridinates to radians:
latNPd = 90 ; % N, degrees
lonNPd = 0 ; % E, degrees
rateNP = 1 ; % unit rate
latNPr = latNPd * conv ; 
lonNPr = lonNPd * conv ;

%Converting our point at 0N (deg) and 0E (deg), with unit radius: 
latPd = 0 ;
lonPd = 0 ;
rateP = 1 ;
latPr = latPd * conv ;
lonPr = lonPd * conv ;

%Converting North Pole spherical cooridinates (radians) to cartesian: 
wxNP = rateNP * cos(lonNPr) * cos(latNPr) ; 
wyNP = rateNP * sin(lonNPr) * cos(latNPr) ;
wzNP = rateNP * sin(latNPr) ; 

%Converting our point spherical cooridinates (radians) to cartesian: 
rxP = rateP * cos(lonPr) * cos(latPr) ; 
ryP = rateP * sin(lonPr) * cos(latPr) ;
rzP = rateP * sin(latPr) ; 

%Calculating absolute motion of Haleakala plate (degrees/my)
Pvxr=(wyNP*rzP)-(wzNP*ryP);
Pvyr=((-wxNP)*rzP)+(wzNP*rxP);
Pvzr=(wxNP*ryP)-(wyNP*rxP);

%Scale from degrees/my to km/my (AKA mm/yr)
sPvx=scale*Pvxr ;
sPvy=scale*Pvyr ;
sPvz=scale*Pvzr ;

%Calculating the length of the velocity 
lengthPv=sqrt(sPvx^2 + sPvy^2 + sPvz^2);

%Printing to answer file
fprintf(outputfile,'\nPoint at 0N, 0E Cartesian Velocity:\n');
fprintf(outputfile,'(using North Pole as rotation pole)\n');
fprintf(outputfile,'vx,      vy,     vz length\n');
fprintf(outputfile,'%0.1f  %9.1f  %6.1f  %6.3f\n',sPvx, sPvy, sPvz, lengthPv);

%%%
% Question 1b
% Calculate the velocity of Haleakala in latitude and longitude. Verify the
% length is the same as 1a.
%%%

% Calculating the velocity in latitude and longitude for Haleakala
latHv = ((-sin(latHr))*(cos(lonHr)*sHvx)) - (sin(latHr)*sin(lonHr)*sHvy) + (cos(latHr)*sHvz);
lonHv = ((-sin(lonHr))*sHvx) + (cos(lonHr)*sHvy);

%Calculating the velocity in latitude and lognitude for our point
latPv = ((-sin(latPr))*(cos(lonPr)*sPvx)) - (sin(latPr)*sin(lonPr)*sPvy) + (cos(latPr)*sPvz);
lonPv = ((-sin(lonPr))*sPvx) + (cos(lonPr)*sPvy);

fprintf(outputfile,'\nQuestion 1b:\n');
fprintf(outputfile,'\nHaleakala Velocity in Lat/Lon:\n');
fprintf(outputfile,'Lat     Lon \n');
fprintf(outputfile,'%0.1f  %4.1f\n',latHv, lonHv);

fprintf(outputfile,'\nOur Point Velocity in Lat/Lon:\n');
fprintf(outputfile,'Lat     Lon \n');
fprintf(outputfile,'%0.1f  %4.1f\n',latPv, lonPv);

%Verifying that the lengths of velocity match
fprintf(outputfile,'\nUsing the UNAVCO Plate Motion Calculator, the\n');
fprintf(outputfile,'speed of Haleakala was verified to be 103 mm/yr.\n');


%%%
% Question 1c
% Calculate the azimuth of the motion vector
%%%

azH = atan2d(lonHv,latHv);
cwazH = 360 + azH;

cwazP = atan2d(lonPv,latPv);

fprintf(outputfile,'\nQuestion 1c:\n');
fprintf(outputfile,'\nHaleakala Azimuth in degrees\n');
fprintf(outputfile,'%0.1f\n',cwazH);
fprintf(outputfile,'\nOur Point Azimuth in degrees\n');
fprintf(outputfile,'%0.1f\n',cwazP);

%%%
% Question 2a
% Calculate a predicted age for the bend in the seamount chain.
%%%

% I already have the latPA_HSr and lonPA_HSr for Haleakala, so now I just
% need to calculate the position of the bend in radians


% I'm using the coordinates of the Daikakuji Seamount as an approximate
% position for the bend
latBd = 32.1 ;
lonBd = 172.3 ;
latBr = latBd * conv ;
lonBr = lonBd * conv ;

%Calculating the distance between Haleakala and the bend:
distHBd = acos((sin(latPA_HSd)*sin(latBd)) + (cos(latPA_HSd)*cos(latBd)*cos(lonBd-lonPA_HSd)))*radius ;
distHBr = acos((sin(latPA_HSr)*sin(latBr)) + (cos(latPA_HSr)*cos(latBr)*cos(lonBr-lonPA_HSr)))*radius ;
% Converting to millimeters
distHBmmd = distHBd * 1000000 ;
distHBmmr = distHBr * 1000000 ;

% Finding the age from absolute velocity
aged = distHBmmd / lengthHv
ager = distHBmmr / lengthHv

%%%
% Question 2b
% (not a matlab question, answer in sheet)
%%%

%%%
% Question 3
% Calculate the vector of Lucknow, India.
%%%

% Latitude
latINEUd = latTSd(27); %degrees
latINEUr = latTSr(27); %radians
% Longitude
lonINEUd = lonTSd(27); %degrees
lonINEUr = lonTSr(27); %radians
% Rate
rateINEU = rateTS(27); %degree/my

% Lucknow location in spherical coordinates (in degrees)
latLd = 26.8 ;
lonLd = 80.9 ;
rateL = 1 ; % using  unit radius (rate information from pole rotation)


% Print INEU and L to HW2.out:
fprintf(outputfile,'\nThe relative plate velocity of the IN-EU Pole: \n');
fprintf(outputfile,'Lat(N), Lon(E),   Rate (deg/my):\n');
fprintf(outputfile,'%1.1f  %5.1f  %9.3f\n',latINEUd, lonINEUd, rateINEU);
fprintf(outputfile,'\nThe location of Lucknow: \n');
fprintf(outputfile,'Lat(N), Lon(E),   Rate (unit rate):\n');
fprintf(outputfile,'%1.1f  %8.1f  %7.3f\n',latLd, lonLd, rateL);

% Calculating the cartesian velocity of Lucknow

% Converting the Lucknow point spherical coordinates to radians:
latLr = latLd * conv ; 
lonLr = lonLd * conv ;

%Converting IN-EU rotation pole coordinates from spherical (radians) into
%cartesian:
wxINEU = rateINEU * cos(lonINEUr) * cos(latINEUr) ; 
wyINEU = rateINEU * sin(lonINEUr) * cos(latINEUr) ;
wzINEU = rateINEU * sin(latINEUr) ; 

%Converting Lucknow point spherical cooridinates (radians) to cartesian: 
wxL = rateL * cos(lonLr) * cos(latLr) ; 
wyL = rateL * sin(lonLr) * cos(latLr) ;
wzL = rateL * sin(latLr) ; 

%Calculating relative motion of Lucknow (degrees/my)
Lvx=(wyINEU*wzL)-(wzINEU*wyL);
Lvy=((-wxINEU)*wzL)+(wzINEU*wxL);
Lvz=(wxINEU*wyL)-(wyINEU*wxL);

%Scale from degrees/my to km/my (AKA mm/yr)
sLvx=scale*Lvx ;
sLvy=scale*Lvy ;
sLvz=scale*Lvz ;

% Calculating the velocity in latitude and longitude for Haleakala
latLv = ((-sin(latLr))*(cos(lonLr)*sLvx)) - (sin(latLr)*sin(lonLr)*sLvy) + (cos(latLr)*sLvz);
lonLv = ((-sin(lonLr))*sLvx) + (cos(lonLr)*sLvy);

fprintf(outputfile,'\nQuestion 3:\n');
fprintf(outputfile,'\nLucknow Velocity in Lat/Lon:\n');
fprintf(outputfile,'Lat     Lon \n');
fprintf(outputfile,'%0.1f  %4.1f\n',latLv, lonLv);

% Calculate the azimuth of the motion vector

azL = atan2d(lonLv,latLv)


fprintf(outputfile,'\nLucknow Azimuth in degrees\n');
fprintf(outputfile,'%0.1f\n',azL);
