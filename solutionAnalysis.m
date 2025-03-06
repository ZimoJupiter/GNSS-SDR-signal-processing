load('navSolCT_1ms_Urban.mat')
urban = navSolutionsCT;
load('navSolCT_1ms_Opensky.mat')
openSky = navSolutionsCT;

urbanGT = [22.3198722, 114.209101777778];
urbanLat = urban.usrPosLLH(:,1);
urbanLon = urban.usrPosLLH(:,2);
urbanLocalTime = urban.localTime(:);

openSkyGT = [22.328444770087565, 114.1713630049711];
openSkyLat = openSky.usrPosLLH(:,1);
openSkyLon = openSky.usrPosLLH(:,2);
openSkyLocalTime = openSky.localTime(:);

figure
geoscatter(urbanLat,urbanLon,2,'r')
hold on
geoscatter(urbanGT(1),urbanGT(2),30,'c','filled')
geoscatter(openSkyLat,openSkyLon,2,'b')
geoscatter(openSkyGT(1),openSkyGT(2),30,'g','filled')
title('Positioning Results')
legend('Urban','Urban GT','Open Sky','Open Sky GT')
hold off

%% Urban
% Lat & Lon
figure
plot(urbanLocalTime,urbanLat)
title('Urban Latitude')
figure
plot(urbanLocalTime,urbanLon)
title('Urban Longitude')

% Ht
figure
plot(urbanLocalTime,urban.usrPosLLH(:,3))
title('Urban Height')

% Vel
figure
plot(urbanLocalTime,urban.usrVelENU)
title('Urban ENU Velocity')
legend('E','N','U')

%% Open Sky
% Lat & Lon
figure
plot(openSkyLocalTime,openSkyLat)
title('Open Sky Latitude')
figure
plot(openSkyLocalTime,openSkyLon)
title('Open Sky Longitude')

% Ht
figure
plot(openSkyLocalTime,openSky.usrPosLLH(:,3))
title('Open Sky Height')

% Vel
figure
plot(openSkyLocalTime,openSky.usrVelENU)
title('Open Sky ENU Velocity')
legend('E','N','U')
