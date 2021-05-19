function [globalForces contactAreas globalPressures ForcesPerSensor AreasPerSensor PressPerSensor] = do_Calculations(Data)
%DOCALCULATIONS Calculates global forces based on tactile data from JOIA data.
%   [GLOBALFORCES CONTACTAREAS PRESSPERSENSOR FORCESPERSENSOR AREASPERSENSOR] = ...
%   DOCALCULATIONS( DATA ) Calculates forces, pressures and contact areas, 
%   based on tactile sensors data from JOIA data previously read by READJOIADATA.m.
%  (outputs are for each of the 4 sensors and the total from all sensors).
%
%   INPUTS:
%   Data: Data that has been read by the READJOIADATA.M function.
% 
%   OUTPUTS:
%   globalForces:    Global forces (in kN) computed from tactile data.
%   contactAreas:    Contact areas  (in m^2) computed from tactile data.
%   globalPressures: Global pressures (in MPa) computed from tactile data.
%   ForcesPerSensor: Forces (in kN) computed from tactile data, for each sensor (from left to right).
%   AreasPerSensor:  Contact areas  (in m^2) computed from tactile data, for each sensor (from left to right).
%   PressPerSensor:  Pressures (in MPa) computed from tactile data, for each sensor (from left to right).
% 

% 
% 
%   Martin Richard
%   C-CORE, St. John's, NL CANADA


%%
% --
% SensorArea = Data.TactileSensor.ASFfile.ROW_SPACING*Data.TactileSensor.ASFfile.COL_SPACING*1e-6; % convert from mm^2 to m^2
SensorArea = Data.TactileSensor.ASFfile.SENSEL_AREA*1e-6; % convert from mm^2 to m^2
SatPressure = Data.TactileSensor.ASFfile.SATURATION_PRESSURE; % MPa
colPerSensor = Data.TactileSensor.ASFfile.COLS/4;

numFrames = size(Data.TactileSensor.ASFfile.Pressures,3);
globalForces    = nan(numFrames,1);
contactAreas    = nan(numFrames,1);
globalPressures = nan(numFrames,1);
ForcesPerSensor	= nan(numFrames,4);
AreasPerSensor  = nan(numFrames,4);
PressPerSensor  = nan(numFrames,4);

didWarn = 0;
for i=1:numFrames
    currentPress = Data.TactileSensor.ASFfile.Pressures(:,:,i);
    if any( abs(currentPress(:)-SatPressure)<0.01 ) && ~didWarn
%         fprintf('------ WARNING: saturation pressure attained for frame %0.0f ------\n',i)
%         fprintf('------ WARNING: Saturation pressure attained for at least one frame ------\n')
        didWarn = 1;
    end
    % NOTE: Data.TactileSensor.ASFfile.Pressures is in MPa - will *1000 to use kPa in order to obtain kN
    localForces     	= currentPress*1000 .* SensorArea;
    globalForces(i,1)   = sum(localForces(:));
    contactAreas(i,1)   = sum(currentPress(:)>0) * SensorArea;
    for iSensor=1:4
        ind = [(iSensor-1)*colPerSensor+1  (iSensor-1)*colPerSensor+colPerSensor];
        currentSensorPress  = currentPress(:,ind(1):ind(2));
        currentSensorForces = localForces(:,ind(1):ind(2));
%         PressPerSensor(i,iSensor)  = sum(currentSensorPress(:));
        ForcesPerSensor(i,iSensor) = sum(currentSensorForces(:));
        AreasPerSensor(i,iSensor)  = sum(currentSensorPress(:)>0) * SensorArea;
    end
end
% -- calculate global pressures
globalPressures = (globalForces./1000) ./ contactAreas;
% -- calculate mean pressures on each panel
for iSensor=1:4
    PressPerSensor(:,iSensor) = (ForcesPerSensor(:,iSensor)./1000) ./ AreasPerSensor(:,iSensor);
end
