function foundHPZs_data = find_HPZs( X,Y,press,handles )
%FINDHPZS retrieve HPZ and stats in frames. Works with JOIATOOL.
% 
%   Ridwan Hossain
%   Memorial University, St. John's, NL CANADA



%%
% NOTES: need (???) to add following:
%   - # of HPZs and/or peaks at each time (i.e., on a given frame)
%   - # of peaks per HPZs at each time
%   - distance between HPZs
%   - shape characteristics of HPZs
%   - Do we keep some consistent tresholds between close frames (i.e., do
%     we change the thresholds between i and i+1 or keep one unique
%     threshold that could be based on the stats of multiple frames as
%     opposed to only one frame ?
% 
% 


%% OPTIONS
normalizationType = 1;      % =1 to center to centroid,
                            % =2 to set scale so that x = [0 1] and y = [0 1] + will compute the circular view
useExpandedVersion = 0;     % determines the number of points that define a polygon 
                            % and therefore the number of points that are used when 
                            % transforming to polar coordinates
                            % =0, uses number of points as determines by contourc.m (runs faster)
                            % =1, uses additional points on the lines (runs slower)
typeOfCorrelationData = 0;  % determines the type of the 'CorrelationData' cell that stores 
                            % information about position (in polar coordinates) of points 
                            % and their neighbours.
                            % = 0, doesn't store anything (faster of the 3 choices)
                            % = 1, stores theta, distance and next distance (results in a 3D data space)
                            % = 2, stores theta, distance, next theta and next distance (slower of the 3 choices slower as it results in a 4D data space)
makePlots = 0;

%% init.
fields2keep = {};
id_2bReseted = [];
% -- contours and shape description
c_contour                   = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'c_contour';
originalContour             = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'originalContour';
normalizedContour           = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'normalizedContour';
normalizedVirtualPolyCenter = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'normalizedVirtualPolyCenter';
theta_dist_relationship     = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'theta_dist_relationship';
distBetweenPoints           = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'distBetweenPoints';
normalizedContour_exp       = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'normalizedContour_exp';         % more points to define the polygon
theta_dist_relationship_exp = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'theta_dist_relationship_exp';   % more points to define the polygon
distBetweenPoints_exp       = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'distBetweenPoints_exp';         % more points to define the polygon
dataFromOutCircle           = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'dataFromOutCircle';             % theta/distances based on an outside circular point of view
% -- ... related to time-dependant processes
fields2keep{end+1,1} = 'correlationData';   id_2bReseted(end+1,1) = 1;
fields2keep{end+1,1} = 'bins_theta';        id_2bReseted(end+1,1) = 1;
fields2keep{end+1,1} = 'bins_D';            id_2bReseted(end+1,1) = 1;
% -- image processing tools to get shapes and stats
isCalculateHPZmask          = false; % TRUE or FALSE - TRUE takes more time - used to get image processing toolbox calcutions
if isCalculateHPZmask
    HPZmask                 = zeros(400,1000);  %#ok<*UNRCH>
    [xmask ymask]           = meshgrid( linspace(X(1,1),X(1,end),size(HPZmask,2))  ,  linspace(Y(1,1),Y(end,1),size(HPZmask,1)) );
    dxmask                  = xmask(1,2)-xmask(1,1); dymask = ymask(2,1)-ymask(1,1);
    properties                  = {};   
    fields2keep{end+1,1} = 'HPZmask';       id_2bReseted(end+1,1) = 1;
    fields2keep{end+1,1} = 'properties';    id_2bReseted(end+1,1) = 1; % will only output something if isCalculateHPZmask==true
else
    HPZmask                 = [];
end
% -- stats related to shape
boundingBox                 = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'boundingBox';
Perimeter                   = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'Perimeter';
InertialMoments             = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'InertialMoments';
CentroidalPrincipalMoments	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'CentroidalPrincipalMoments';
HPZareas                  	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'HPZareas';
HPZcentroids             	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'HPZcentroids';
num_sensorsInHPZ            = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'num_sensorsInHPZ';
% -- threshold values
if get(handles.saveType,'SelectedObject')==handles.varyingThreshold
    isConstantThreshold = 0;
    error = 0; try mult_std_peaks  = str2double(get(handles.mult_std_peaks,'String')); catch, error = 1; end %#ok<*CTCH>
    if error || ~isnumeric(mult_std_peaks) || isempty(mult_std_peaks) || isnan(mult_std_peaks), mult_std_peaks = 1.0; end
    error = 0; try mult_std_limits = str2double(get(handles.mult_std_limits,'String')); catch, error = 1; end
    if error || ~isnumeric(mult_std_limits) || isempty(mult_std_limits) || isnan(mult_std_limits), mult_std_limits = 2.0; end
    error = 0; try minValue_peaks = str2double(get(handles.minValue_peaks,'String')); catch, error = 1; end
    if error || ~isnumeric(minValue_peaks) || isempty(minValue_peaks) || isnan(minValue_peaks), minValue_peaks = 3.0; end
    fields2keep{end+1,1} = 'mult_std_peaks';    id_2bReseted(end+1,1) = 0;
    fields2keep{end+1,1} = 'mult_std_limits';   id_2bReseted(end+1,1) = 0;
    fields2keep{end+1,1} = 'minValue_peaks';    id_2bReseted(end+1,1) = 0;
else
    isConstantThreshold = 1;
    error = 0; try ConstantThreshold  = str2double(get(handles.edit19,'String')); catch, error = 1; end %#ok<*CTCH>
    if error || ~isnumeric(ConstantThreshold) || isempty(ConstantThreshold) || isnan(ConstantThreshold), ConstantThreshold = 2.0; end
    fields2keep{end+1,1} = 'ConstantThreshold'; id_2bReseted(end+1,1) = 0;
end
treshold4peaks              = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'treshold4peaks'; %#ok<*NASGU>
treshold4limits             = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'treshold4limits';
% -- ... related to pressures and thresholds for determining HPZs
peaks                       = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'peaks';
allPeaks                   	= []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'allPeaks';
num_peaksInHPZ              = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'num_peaksInHPZ';
closestPeak                 = {}; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'closestPeak';

meanHPZpressures           	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'meanHPZpressures'; % -- defined for all threshold tested
stdHPZpressures           	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'stdHPZpressures'; % -- defined for all threshold tested
maxHPZpressures           	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'maxHPZpressures'; % -- defined for all threshold tested
minHPZpressures           	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'minHPZpressures'; % -- defined for all threshold tested
glob_press_onHPZ            = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'glob_press_onHPZ'; % -- defined for all threshold tested
Force_HPZ                   = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'Force_HPZ'; % IN NEWTON  -- defined for all threshold tested
PressGETthreshold       	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'PressGETthreshold'; % pressures greater or equal than threshold -- defined for all threshold tested
PressLTthreshold            = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'PressLTthreshold';  % pressures lower or equal than threshold -- defined for all threshold tested
MeanPressGETthreshold   	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'MeanPressGETthreshold'; % -- defined for all threshold tested
MeanPressLTthreshold        = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'MeanPressLTthreshold'; % -- defined for all threshold tested
assocForceGETthreshold   	= []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'assocForceGETthreshold'; % -- defined for all threshold tested
assocForceLTthreshold       = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'assocForceLTthreshold'; % -- defined for all threshold tested
areaAssocGETthreshold       = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'areaAssocGETthreshold'; % -- defined for all threshold tested
areaAssocLTthreshold        = []; id_2bReseted(end+1,1) = 1; fields2keep{end+1,1} = 'areaAssocLTthreshold'; % -- defined for all threshold tested

thresholdsSequence          = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'thresholdsSequence'; % -- defined for all threshold tested
areasForGivenThreshold      = {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'areasForGivenThreshold'; % -- defined for all threshold tested
centroidsForGivenThreshold	= {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'centroidsForGivenThreshold'; % -- defined for all threshold tested
thresholdsAreaRelationships = {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'thresholdsAreaRelationships'; % -- defined for all threshold tested
meanHPZpressures_ATR        = {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'meanHPZpressures_ATR'; % -- defined for all threshold tested
meanHPZpressures_ATRtmp     = {};
stdHPZpressures_ATR       	= {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'stdHPZpressures_ATR'; % -- defined for all threshold tested
stdHPZpressures_ATRtmp      = {};
maxHPZpressures_ATR      	= {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'maxHPZpressures_ATR'; % -- defined for all threshold tested
maxHPZpressures_ATRtmp      = {};
minHPZpressures_ATR     	= {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'minHPZpressures_ATR'; % -- defined for all threshold tested
minHPZpressures_ATRtmp      = {};
glob_press_onHPZ_ATR    	= {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'glob_press_onHPZ_ATR'; % -- defined for all threshold tested
Force_HPZ_ATR            	= {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'Force_HPZ_ATR'; % -- defined for all threshold tested
Force_HPZ_ATRtmp            = {};
PressGETthreshold_ATR       = {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'PressGETthreshold_ATR'; % -- defined for all threshold tested
PressLTthreshold_ATR        = {}; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'PressLTthreshold_ATR'; % -- defined for all threshold tested
MeanPressGETthreshold_ATR   = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'MeanPressGETthreshold_ATR'; % -- defined for all threshold tested
MeanPressLTthreshold_ATR    = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'MeanPressLTthreshold_ATR'; % -- defined for all threshold tested
assocForceGETthreshold_ATR  = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'assocForceGETthreshold_ATR'; % -- defined for all threshold tested
assocForceLTthreshold_ATR   = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'assocForceLTthreshold_ATR'; % -- defined for all threshold tested
areaAssocGETthreshold_ATR   = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'areaAssocGETthreshold_ATR'; % -- defined for all threshold tested
areaAssocLTthreshold_ATR    = []; id_2bReseted(end+1,1) = 0; fields2keep{end+1,1} = 'areaAssocLTthreshold_ATR'; % -- defined for all threshold tested

% --


%% probability distributions between points of polygons defining shapes - prepare bins for theta and distance
% -- bins for distance
numDBins = 55;  % number of bins for distance
maxD = 110;     % maximum distance for bins creation
all_D = linspace(0,maxD,numDBins+1);        % all distances   
bins_D = [all_D(1:end-1)' all_D(2:end)'];	% bins ( [lower & upper[ limits)
% -- bins for theta
%	(the larger the window - i.e., smaller is the number of bins - the more randomness it adds to generated shapes)
numThetaBins = 120;	% number of bins for theta
all_thetas = linspace(-pi,pi,numThetaBins+1);     % all thetas 
bins_theta = [all_thetas(1:end-1)' all_thetas(2:end)'];	% bins ( [lower & upper[ limits)
bins_theta(end,2) = bins_theta(end,2)+pi/360; % needs to be slightly >pi because the condition used lower is <bins_theta(:,2) and as this is the last one we should use <=bins_theta(end,2), so this bypass the potential issue
% -- prepare 'correlationData' cell into which the data will be stored later
if typeOfCorrelationData==0
    correlationData = {};
elseif typeOfCorrelationData==1
    %	'correlationData' contains theta->distance->next_distance
    correlationData     = cell( size(bins_theta,1),1 );
    blank_distCell      = cell( size(bins_D,1),1 );
    null_distMatrix     = zeros( size(bins_D,1),1 );
    for iCorr = 1:numel(correlationData)
        correlationData{iCorr} = blank_distCell;
        for jCorr = 1:size(bins_D,1)
            correlationData{iCorr}{jCorr,1} = null_distMatrix;
        end
    end
elseif typeOfCorrelationData==2
    % -- prepare 'correlationData2' cell into which the data will be stored later
    %	'correlationData' contains theta->distance->next_theta->next_distance
    blank_thetaCell   	= cell( size(bins_theta,1),1 );
    correlationData 	= blank_thetaCell;
    blank_distCell      = cell( size(bins_D,1),1 );
    null_distMatrix     = zeros( size(bins_D,1),1 );
    for iCorr = 1:numel(correlationData)
        correlationData{iCorr} = blank_distCell;
        for jCorr = 1:size(bins_D,1)
            correlationData{iCorr}{jCorr,1} = blank_thetaCell;
            for kCorr = 1:size(blank_thetaCell,1)
                correlationData{iCorr}{jCorr,1}{kCorr,1} = null_distMatrix;
            end
        end
    end
end
% --


%% FIND HPZ & LOOP THROUGH EACH CALCULATED CONTOUR
% -- find all peaks and positive pressures
[~,imax,~,~] = extrema2(press); % indices for all peaks
allPeaks = [X(imax) Y(imax) press(imax)]; % get remaining peaks coordinates
ind = press>0;
PositivePressures = press(ind); % get all positive pressures (rid of zeros for stats calculations)
if isConstantThreshold
    treshold4limits = ConstantThreshold;
    treshold4peaks  = ConstantThreshold;
%     allThresholdsValues2Test = [1:0.25:7];
    allThresholdsValues2Test = [];
else
    treshold4limits = max( [ mean(PositivePressures(:)) + mult_std_peaks *std(PositivePressures(:))   minValue_peaks  ] ); % in MPa
    treshold4peaks  = max( [ mean(PositivePressures(:)) + mult_std_limits*std(PositivePressures(:))   treshold4limits ] ); % in MPa
    allThresholdsValues2Test = [];
end


% -- DO THE SEQUENCE OF THRESHOLD IF ASKED
if ~isempty(allThresholdsValues2Test)
    
    % --
    % PACK all variables in the workspace as of now
    AllVarIN = whos;
    for iIN = 1:numel(AllVarIN)
        if ~ismember(AllVarIN(iIN).name,{'AllVar','AllVarIN','AllVarOUT'})
            AllVarIN(iIN).Data = eval(AllVarIN(iIN).name);
        end
    end
    % --
    for treshold4limits=allThresholdsValues2Test
        treshold4peaks = treshold4limits;
        
        % -- change the values of some important variables into the AllVarIN structure
        for iIN = 1:numel(AllVarIN)
            switch AllVarIN(iIN).name
                case {'treshold4limits','treshold4peaks','thresholdsSequence','areasForGivenThreshold','centroidsForGivenThreshold',...
                      'assocForceGETthreshold_ATR','assocForceLTthreshold_ATR','meanHPZpressures_ATR','stdHPZpressures_ATR',...
                      'maxHPZpressures_ATR','minHPZpressures_ATR','PressGETthreshold_ATR','PressLTthreshold_ATR',...
                      'MeanPressGETthreshold_ATR','MeanPressLTthreshold_ATR','areaAssocGETthreshold_ATR','areaAssocLTthreshold_ATR',...
                      'Force_HPZ_ATR','meanHPZpressures_ATRtmp','stdHPZpressures_ATRtmp','maxHPZpressures_ATRtmp',...
                      'minHPZpressures_ATRtmp','Force_HPZ_ATRtmp'};
                    AllVarIN(iIN).Data = eval(AllVarIN(iIN).name);
            end
        end
        
        % call the routine
        AllVarOUT = findHPZ_and_getHPZprops(AllVarIN);
        % UNPACK the variables in the workspace
        for iOUT = 1:numel(AllVarOUT)
            if ~ismember(AllVarOUT(iOUT).name,{'AllVar','AllVarIN','AllVarOUT'})
                eval( [AllVarOUT(iOUT).name,' = AllVarOUT(iOUT).Data;'] )
            end
        end
        % --
        
        % -- place data in here ...
        thresholdsSequence(end+1,1)         = treshold4limits;
        areasForGivenThreshold{end+1,1}     = HPZareas;
        centroidsForGivenThreshold{end+1,1} = HPZcentroids;
        meanHPZpressures_ATRtmp{end+1,1}  	= meanHPZpressures;
        stdHPZpressures_ATRtmp{end+1,1}   	= stdHPZpressures;
        maxHPZpressures_ATRtmp{end+1,1}   	= maxHPZpressures;
        minHPZpressures_ATRtmp{end+1,1}   	= minHPZpressures;
        Force_HPZ_ATRtmp{end+1,1}         	= Force_HPZ;
%         PressGETthreshold_ATR{end+1,1}     	= PressGETthreshold;
%         PressLTthreshold_ATR{end+1,1}       = PressLTthreshold;
        MeanPressGETthreshold_ATR(end+1,1) 	= MeanPressGETthreshold;
        MeanPressLTthreshold_ATR(end+1,1)  	= MeanPressLTthreshold;
        assocForceGETthreshold_ATR(end+1,1) = assocForceGETthreshold;
        assocForceLTthreshold_ATR(end+1,1)	= assocForceLTthreshold;
        areaAssocGETthreshold_ATR(end+1,1)	= areaAssocGETthreshold;
        areaAssocLTthreshold_ATR(end+1,1)  	= areaAssocLTthreshold;
        % --

    end
    
    % -- need to track the centroids to come up with the curves for each HPZ
    pos2track = [];
    for iStep = 1:numel(thresholdsSequence)
      	centroid = centroidsForGivenThreshold{iStep};
        if ~isempty(centroid)
            pos2track = [pos2track ; centroid  ones(size(centroid,1),1)*iStep ];
        end
    end
    % -- call tracking routine
    try
        maxDisplacement = 30;
        resTracking = track(pos2track,maxDisplacement); % output is: [ (x) (y) (t) (id) ]
    catch
        maxDisplacement = 20;
        resTracking = track(pos2track,maxDisplacement); % output is: [ (x) (y) (t) (id) ]
    end
    % -- construct curves
    individualHPZ = unique(resTracking(:,end));
    % -- find out which particles are 'alone', i.e., appear only in one time step
    HPZ_alone = [];
    for i=1:numel(individualHPZ)
        if sum(individualHPZ(i)==resTracking(:,end))==1
            HPZ_alone(end+1,1) = individualHPZ(i);
        end
    end
    % --
    for i=1:numel(individualHPZ) % for all individual HPZs
        if ~ismember( individualHPZ(i),HPZ_alone )
            ind                             = resTracking(:,end)==individualHPZ(i);
            centroidsTimeTrace              = resTracking(ind,1:end-2);
            PositionInThresholdsSequence  	= resTracking(ind,end-1);
            % -- find out which data is associated with it
            currentTimesTraces = [];
            iCount = 1;
            for iPos = PositionInThresholdsSequence(1):PositionInThresholdsSequence(end)
                centroidsOfNeededPos = centroidsForGivenThreshold{iPos};
                dist = sqrt((centroidsOfNeededPos(:,1)-centroidsTimeTrace(iCount,1)).^2  +  (centroidsOfNeededPos(:,2)-centroidsTimeTrace(iCount,2)).^2);
                ind_dist = find(dist==0, 1); % this will be the ID specific to the curent frame
                if ~isempty(ind_dist)
                    current_thresh              = thresholdsSequence(iPos);
                    current_area                = areasForGivenThreshold{iPos}(ind_dist);
                    meanHPZpressures_tmp     	= meanHPZpressures_ATRtmp{iPos}(ind_dist);
                    stdHPZpressures_tmp     	= stdHPZpressures_ATRtmp{iPos}(ind_dist);
                    maxHPZpressures_tmp     	= maxHPZpressures_ATRtmp{iPos}(ind_dist);
                    minHPZpressures_tmp     	= minHPZpressures_ATRtmp{iPos}(ind_dist);
                    Force_HPZ_tmp               = Force_HPZ_ATRtmp{iPos}(ind_dist);
                    currentTimesTraces(end+1,:) = [ current_thresh  current_area ...
                                                    meanHPZpressures_tmp  stdHPZpressures_tmp  maxHPZpressures_tmp  minHPZpressures_tmp ...
                                                    Force_HPZ_tmp
                                                  ];
                end
                iCount = iCount + 1;
            end
            thresholdsAreaRelationships{end+1,1}    = currentTimesTraces(:,1:2);
            meanHPZpressures_ATR{end+1,1}           = currentTimesTraces(:,3);
         	stdHPZpressures_ATR{end+1,1}            = currentTimesTraces(:,4);
          	maxHPZpressures_ATR{end+1,1}            = currentTimesTraces(:,5);
        	minHPZpressures_ATR{end+1,1}            = currentTimesTraces(:,6);
        	Force_HPZ_ATR{end+1,1}                  = currentTimesTraces(:,7);
        end
    end
    
    % -- RESET VARIABLES THAT WERE SETUP AND WEREN'T SUPPOSED TO
    for i=1:numel(id_2bReseted)
        if id_2bReseted(i)
            if eval( ['iscell(',fields2keep{i},')'] )
                eval( [fields2keep{i},' = {};'] )
            else
                eval( [fields2keep{i},' = [];'] )
            end
        end
    end
    % -- 
    
end


% -- DO THE ASKED THRESHOLD IN THE END ONLY --
% PACK all variables in the workspace as of now
AllVarIN = whos;
for iIN = 1:numel(AllVarIN)
    if ~ismember(AllVarIN(iIN).name,{'AllVar','AllVarIN','AllVarOUT'})
        AllVarIN(iIN).Data = eval(AllVarIN(iIN).name);
    end
end
% call the routine
AllVarOUT = findHPZ_and_getHPZprops(AllVarIN);
% UNPACK the variables in the workspace BUT NOT THOSE DONE BEFORE IN THE PREVIOUS BLOCK
for iOUT = 1:numel(AllVarOUT)
    % --
    doUnpack = true;
    % --
    if ismember(AllVarOUT(iOUT).name,{'fields2keep','id_2bReseted','AllVar','AllVarIN','AllVarOUT'})
        doUnpack = false;
    elseif ismember(AllVarOUT(iOUT).name,fields2keep)
        ind = strfind(fields2keep,AllVarOUT(iOUT).name);
        for iTmp = 1:numel(ind)
            if ~isempty( ind{iTmp} )
                ind = iTmp;
                break;
            end
        end
        if ~id_2bReseted(ind)
            doUnpack = false;
        end
    end
    if doUnpack
        eval( [AllVarOUT(iOUT).name,' = AllVarOUT(iOUT).Data;'] )
    end
end
% --



%% keep desired data and output
% --
for i = 1:numel(fields2keep)
    eval( ['foundHPZs_data.(fields2keep{i}) = ',fields2keep{i},';'] )
end
% --


%%
function AllVarOUT = findHPZ_and_getHPZprops(AllVar)
% retrieve the variables that were in the workspace
for iIN = 1:numel(AllVar)
    if ~ismember(AllVar(iIN).name,{'AllVar','AllVarIN','AllVarOUT'})
        eval( [AllVar(iIN).name,' = AllVar(iIN).Data;'] )
    end
end
% --

% FINDS CONTOURS OF HPZ BASED ON THRESHOLD
% --
flags2erase = press(imax)<treshold4peaks; % flags to erase all peaks that are lower than the peak threshold
imax(flags2erase) = []; % ... erase
peaks = [X(imax) Y(imax) press(imax)]; % get remaining peaks coordinates
% -- calculte HPZs contours, based on calculated thresholds
c_contour_temp = contourc(X(1,:),Y(:,1),press,[treshold4limits treshold4limits]); % NOTE: Use CONTOURC(X, Y, Z, [v, v]) to compute a single contour at the level v.
% --

% CALCULATES PROPERTIES OF EACH FOUND HPZ
% --
AllInd_HPZ = false(size(X));
unitArea = (5.19938/1000)^2;
while ~isempty(c_contour_temp)
    % -- get current contour
    num_elem = c_contour_temp(2,1);
    currentContour = c_contour_temp(:,2:num_elem+1)';
    if ~(currentContour(end,1)==currentContour(1,1) && currentContour(end,2)==currentContour(1,2))
        currentContour = [currentContour ; currentContour(1,:)]; % to close the polygon in that case
    end
    % -- check if the current HPZ is at least bigger than the minimum allowed area for a HPZ - if not, reject (i.e., do not consider it a HPZ)
    TempArea = polyarea(currentContour(:,1),currentContour(:,2));
    if TempArea>5 % in mm^2
        % -- check if there is at least one peak inside the contour - if not, reject (i.e., do not consider it a HPZ)
        in = inpolygon(peaks(:,1),peaks(:,2),currentContour(:,1),currentContour(:,2));
        if ~any(in) % no peaks inside the contour, reject
            currentContour = [];
        else % at least one peak, keep as a HPZ
            % -- get the position of the first peak inside the countour
            firstPeakPos = peaks(find(in,1,'first'),1:2);
            % -- compute stats about the HPZ shape and position
            HPZareas(end+1,1)                   = TempArea; %#ok<*AGROW>
            [geom,IM,CPM]                       = polygeom( currentContour(:,1),currentContour(:,2) );
            HPZcentroids(end+1,:)               = geom(2:3);
            Perimeter(end+1,:)                  = geom(4);
            InertialMoments(end+1,:)            = IM;  % [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ], u,v are centroidal axes parallel to x,y axes.
            CentroidalPrincipalMoments(end+1,:)	= CPM; % [ I1     ang1   I2     ang2   J ], I1,I2 are centroidal principal moments about axes at angles ang1,ang2. ang1 and ang2 are in radians. J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv
            boundingBox(end+1,:)                = [ min(currentContour)  max(currentContour) ]; % [minX minY maxX maxY]
                in = inpolygon(X,Y,currentContour(:,1),currentContour(:,2)); % find indices in the sensor layout that are identified as being the HPZ
            num_sensorsInHPZ(end+1,1)           = sum(in(:));
            meanHPZpressures(end+1,1)           = mean(press(in));
            stdHPZpressures(end+1,1)            = std(press(in));
            maxHPZpressures(end+1,1)            = max(press(in));
            minHPZpressures(end+1,1)            = min(press(in));
            glob_press_onHPZ(end+1,1)           = sum( press(in).*1e6.*unitArea ) / (num_sensorsInHPZ(end).*1e6.*unitArea);
            Force_HPZ(end+1,1)                  = sum( press(in).*1e6.*unitArea ); % NEWTON
                in2 = inpolygon(allPeaks(:,1),allPeaks(:,2),currentContour(:,1),currentContour(:,2));
            num_peaksInHPZ(end+1,1)           = sum(in2(:));
            currentpeakind = find(in2);
            AllInd_HPZ(in) = 1;
            closestPeakTemp = nan(numel(currentpeakind),3);
            if numel(currentpeakind)==1
                closestPeakTemp = NaN;
            else
                for ipeaks = 1:numel(currentpeakind)
                    closestDist = sqrt((allPeaks(currentpeakind(ipeaks),1)-allPeaks(currentpeakind,1)).^2 + (allPeaks(currentpeakind(ipeaks),2)-allPeaks(currentpeakind,2)).^2);
                    closestDist(ipeaks) = [];
                    closestPeakTemp(ipeaks,:) = [ allPeaks(ipeaks,1) allPeaks(ipeaks,2) min(closestDist) ];
                end
            end
            closestPeak{end+1,1} = closestPeakTemp;            
            % -- find the 'center' of the polygon (although technically this does not exist). We need a reference point inside the polygon
            if inpolygon( HPZcentroids(end,1),HPZcentroids(end,2) ,currentContour(:,1),currentContour(:,2) );
                % if centroid inside polygon, use this ...
                virtualPolyCenter = HPZcentroids(end,:);
            else
                % ... otherwise, use the first peak position
                virtualPolyCenter = firstPeakPos;
%                 virtualPolyCenter = my_irregularPolyCenter(currentContour(:,1),currentContour(:,2)); pause; close;
            end  
            % -- keep original shape + do normalization of the HPZ shape
            originalContour{end+1,1} = currentContour;
            if normalizationType==1 % translation so that centroid is at (0,0)
                normalizedContourTemp = [currentContour(:,1)-HPZcentroids(end,1)  currentContour(:,2)-HPZcentroids(end,2)];
                normalizedVirtualPolyCenter(end+1,:) = [0 0];
            elseif normalizationType==2 % translation and scaling so that HPZ is from x = [0 1] and y = [0 1]
                normalizedContourTemp = [currentContour(:,1)-min(currentContour(:,1))  currentContour(:,2)-min(currentContour(:,2))];
                multx = max( normalizedContourTemp(:,1) ) - min( normalizedContourTemp(:,1) ); if multx==0, multx=1; end
                multy = max( normalizedContourTemp(:,2) ) - min( normalizedContourTemp(:,2) ); if multy==0, multy=1; end
                normalizedContourTemp(:,1) = normalizedContourTemp(:,1) ./ multx;
                normalizedContourTemp(:,2) = normalizedContourTemp(:,2) ./ multy;          
                normalizedVirtualPolyCenter(end+1,:) = [ (virtualPolyCenter(1)-min(currentContour(:,1)))/multx  (virtualPolyCenter(2)-min(currentContour(:,2)))/multy ];
            end
            normalizedContour{end+1,1} = normalizedContourTemp;
            % -- description of the shape in terms of polar coordinates (theta and distance)
            [theta,polar_dist] = cart2pol(normalizedContourTemp(:,1),normalizedContourTemp(:,2));
            theta_dist_relationship{end+1,1} = [theta polar_dist];
            distBetweenPoints{end+1,1} = [NaN ; sqrt((normalizedContourTemp(2:end,1)-normalizedContourTemp(1:end-1,1)).^2  +  (normalizedContourTemp(2:end,2)-normalizedContourTemp(1:end-1,2)).^2)];
            % -- optional plots
            if makePlots
                fig0 = figure('Position',[-1612 428 560 420]); hold on; plot(normalizedContourTemp(:,1),normalizedContourTemp(:,2),'o-');
                fig1 = figure('Position',[-1612 428 560 420]); hold on; plot(normalizedContourTemp(:,1),normalizedContourTemp(:,2),'o-'); dx = get(gca,'XLim'); dx = dx(2)-dx(1); dy = get(gca,'YLim'); dy = dy(2)-dy(1); plot(-0.2,-0.2,'or','MarkerFaceColor','r'); for ii=1:size(normalizedContourTemp,1)-1, text(normalizedContourTemp(ii,1)+dx*0.02,normalizedContourTemp(ii,2)*1.03,num2str(ii)); plot([-0.2 normalizedContourTemp(ii,1)],[-0.2 normalizedContourTemp(ii,2)],'-r'); end
                fig2 = figure('Position',[-1000 428 560 420]); plot(theta,polar_dist,'o'); xlabel('\theta'); ylabel('dist.'); dx = get(gca,'XLim');  dx = dx(2)-dx(1); dy = get(gca,'YLim'); dy = dy(2)-dy(1);   for ii=1:size(normalizedContourTemp,1)-1, text(theta(ii)+dx*0.02,polar_dist(ii)+dy*0.03,num2str(ii)); end
                pause; 
                close([fig0 fig1 fig2]);
            end
            % -- use more points to describe the polygon - will add points
            if useExpandedVersion
                normalizedContourTemp_2 = addPointsToPolygon( normalizedContourTemp );               
                normalizedContour_exp{end+1,1} = normalizedContourTemp_2;
                % -- optional plots
                if makePlots
                    fig0 = figure('Position',[-1612 428 560 420]); hold on;
                    plot(normalizedContourTemp(:,1),normalizedContourTemp(:,2),'o-');
                    plot(normalizedContourTemp_2(:,1),normalizedContourTemp_2(:,2),'xr')
                    pause; close(fig0);
                end
                % -- description of the 'expanded' shape in terms of polar coordinates (theta and distance)
                [theta,polar_dist] = cart2pol(normalizedContourTemp_2(:,1),normalizedContourTemp_2(:,2));
                    indErase = isnan(theta)|isnan(polar_dist);
                theta( indErase ) = [];  polar_dist( indErase ) = [];
                [theta inSorted] = sort(theta); polar_dist = polar_dist(inSorted); % sort theta - really important if we want the shape reconstruction process to work
                theta_dist_relationship_exp{end+1,1} = [theta polar_dist];
                distBetweenPoints_exp{end+1,1} = [NaN ; sqrt((normalizedContourTemp_2(2:end,1)-normalizedContourTemp_2(1:end-1,1)).^2  +  (normalizedContourTemp_2(2:end,2)-normalizedContourTemp_2(1:end-1,2)).^2)];
                [theta indSorted] = sort(theta); polar_dist = polar_dist(indSorted);
                if typeOfCorrelationData>0
                    % -- place theta and distances into bins to get statistics that are useable later to reproduce shapes
                    for iBins = 1:numel(polar_dist)-1
                        D_binNo         = find( polar_dist(iBins)>=bins_D(:,1)   &  polar_dist(iBins)<bins_D(:,2)   );  % current dist. bin #
                        theta_binNo     = find( theta(iBins)>=bins_theta(:,1)    &  theta(iBins)<bins_theta(:,2)    );  % current theta bin #
                        nextD_binNo     = find( polar_dist(iBins+1)>=bins_D(:,1) &  polar_dist(iBins+1)<bins_D(:,2) );  % next pint dist. bin #
                        if typeOfCorrelationData==1
                            correlationData{theta_binNo}{D_binNo}(nextD_binNo) = correlationData{theta_binNo}{D_binNo}(nextD_binNo) + 1; % update the count to "count=count+1"
                        elseif typeOfCorrelationData==2
                            nexttheta_binNo = find( theta(iBins+1)>=bins_theta(:,1)  &  theta(iBins+1)<bins_theta(:,2)  );  % next theta bin #
                            correlationData{theta_binNo}{D_binNo}{nexttheta_binNo}(nextD_binNo) = correlationData{theta_binNo}{D_binNo}{nexttheta_binNo}(nextD_binNo) + 1; % update the count to "count=count+1"
                        end
                    end
                end
            end
            % -- use a circle around the normalized shape to calculte stats.
            if normalizationType==2
                allDist = getStatsFromOusideCircle( normalizedContourTemp );
                dataFromOutCircle{end+1,1} = allDist;
                dataFromOutCircle{end,2} = [multx multy];
            end
            % -- use image processing tools to calculte statistics, if wanted (this takes longer to run - obsolete)
            if isCalculateHPZmask
                [HPZmask propertiesTemp] = CalculateHPZmask( currentContour,HPZmask,dxmask,xmask,ymask,HPZareas );
                properties{end+1,1} = propertiesTemp;
            end
            % --

            % -- update contours
            currentContour = [currentContour; [nan nan]];
        end
    else
        currentContour = [];
    end
    % store data inside main c_contour matrix
    c_contour = [ c_contour ; currentContour ];
    c_contour_temp(:,1:num_elem+1) = [];   
end
% --

% store other data related to this particular frame
PressGETthreshold       = press(AllInd_HPZ);
MeanPressGETthreshold	= mean(press(AllInd_HPZ));
assocForceGETthreshold	= sum(press(AllInd_HPZ).*1e6.*unitArea);
areaAssocGETthreshold   = sum(AllInd_HPZ(:)).*unitArea;

PressLTthreshold        = press(~AllInd_HPZ);
MeanPressLTthreshold	= mean(press(~AllInd_HPZ));
assocForceLTthreshold	= sum(press(~AllInd_HPZ).*1e6.*unitArea);
areaAssocLTthreshold 	= sum(~AllInd_HPZ(:)).*unitArea;

% --
if isempty(num_sensorsInHPZ) && (~isempty(treshold4peaks)||~isempty(treshold4limits)||~isempty(peaks))
    treshold4peaks  = [];
    treshold4limits = [];
    peaks           = [];
end

% --
% PACK all variables in the workspace as of now
AllVarOUT = whos;
for iOUT = 1:numel(AllVarOUT)
    if ~ismember(AllVarOUT(iOUT).name,{'AllVar','AllVarIN','AllVarOUT'})
        AllVarOUT(iOUT).Data = eval(AllVarOUT(iOUT).name);
    end
end
% --

%%
function contourWithMorePoints = addPointsToPolygon( originalContour )
delta_xy = 0.1; % in mm
roundingFactor = 1/delta_xy;
contourWithMorePoints(1,:) = originalContour(1,:);
for iPoints = 2:size(originalContour,1)
    % get sign of progression with respect to x (i.e., slope w/r to x)
    sign = double(originalContour(iPoints-1,1)<originalContour(iPoints,1)) * 2 - 1;
    % add points in x at each delta_xy interval (rounded to the closest 1/delta_xy value)
    if sign==1
        nextRoundedValue = ceil(originalContour(iPoints-1,1)*roundingFactor)/roundingFactor;
    elseif sign==-1
        nextRoundedValue = floor(originalContour(iPoints-1,1)*roundingFactor)/roundingFactor;
    end
    xtemp = [ nextRoundedValue:delta_xy*sign:originalContour(iPoints,1)   originalContour(iPoints,1)]';
    % total Dx between the extremities of the line
    xtempDiff = diff( originalContour(iPoints-1:iPoints,1) );
    if any(~xtempDiff)
        sign = double(originalContour(iPoints-1,2)<originalContour(iPoints,2)) * 2 - 1;
            nextRoundedValue = ceil(originalContour(iPoints-1,2)*roundingFactor)/roundingFactor;
        ytemp = [ nextRoundedValue:delta_xy*sign:originalContour(iPoints,2)   originalContour(iPoints,2)]';
        xtemp = interp1( originalContour(iPoints-1:iPoints,2) , originalContour(iPoints-1:iPoints,1) , ytemp );
    else
        ytemp = interp1( originalContour(iPoints-1:iPoints,1) , originalContour(iPoints-1:iPoints,2) , xtemp );
    end
    contourWithMorePoints = [contourWithMorePoints ; [xtemp ytemp]];
end

%%
function allDist = getStatsFromOusideCircle( normalizedContourTemp )
% --
outsideExc = 0.0;
caseWithRectange = 0; % 1 for rectangle, 0 for circle
doplot = 0;
if doplot
    fig3 = figure('Position',[-1612 428 560 420]); hold on; axis equal; grid on;
    plot(normalizedContourTemp(:,1),normalizedContourTemp(:,2),'o-');
    dx = get(gca,'XLim'); dx = dx(2)-dx(1); dy = get(gca,'YLim'); dy = dy(2)-dy(1);
end
centerpt = [0.5 0.5];
% --
if caseWithRectange
    % --- case with a rectangle
    increment = 0.05;
    pts2scan = [];
    temp = [-outsideExc:increment:1+outsideExc]';                       pts2scan = [pts2scan ; temp                                 ones(numel(temp),1)*-outsideExc     ]; %#ok<*NBRAK>
    temp = [-outsideExc+increment:increment:1+outsideExc]';             pts2scan = [pts2scan ; ones(numel(temp),1)*1+outsideExc     temp                                ];
    temp = [1+outsideExc-increment:-increment:-outsideExc]';            pts2scan = [pts2scan ; temp                                 ones(numel(temp),1)*1+outsideExc  	];
    temp = [1+outsideExc-increment:-increment:-outsideExc+increment]';	pts2scan = [pts2scan ; ones(numel(temp),1)*-outsideExc      temp                                ];
    % pts2scan = [  [-0.2:increment:1.2]'   ones(8,1)*-0.2 ;
    %               ones(7,1)*1.2     [0.0:0.2:1.2]' ;
    %               [1.0:-0.2:-0.2]'  ones(7,1)*1.2  ;
    %               ones(6,1)*-0.2    [1.0:-0.2:0]'  ];
    if doplot
        plot([-outsideExc 1+outsideExc 1+outsideExc -outsideExc -outsideExc] , [-outsideExc -outsideExc 1+outsideExc 1+outsideExc -outsideExc],'o-r','MarkerFaceColor','r');
        plot(pts2scan(:,1),pts2scan(:,2),'ok')
    end
    oppositePts = [];
    temp = [1+outsideExc:-increment:-outsideExc]';                      oppositePts = [oppositePts ; temp                               ones(numel(temp),1)*1+outsideExc 	];
    temp = [1+outsideExc-increment:-increment:-outsideExc]';            oppositePts = [oppositePts ; ones(numel(temp),1)*-outsideExc    temp                                ];
    temp = [-outsideExc+increment:increment:1+outsideExc]';             oppositePts = [oppositePts ; temp                               ones(numel(temp),1)*-outsideExc  	];
    temp = [-outsideExc+increment:increment:1+outsideExc-increment]';  	oppositePts = [oppositePts ; ones(numel(temp),1)*1+outsideExc	temp                                ];
    % oppositePts = [  [1.2:-0.2:-0.2]'   ones(8,1)*1.2    ;
    %                  ones(7,1)*-0.2     [1.0:-0.2:-0.2]' ;
    %                  [0.0:0.2:1.2]'     ones(7,1)*-0.2  ;
    %                  ones(6,1)*1.2      [0.0:0.2:1]'  ];
    if doplot
        set(gca,'XLim',[-outsideExc-0.01 1+outsideExc+0.01],'YLim',[-outsideExc-0.01 1+outsideExc+0.01])
    end
else
    % --- case with a circle
    radius = sqrt( (1+outsideExc-centerpt(2))^2 + (0-outsideExc-centerpt(1))^2 );
    theta = [0:1:359]';
    pts2scan = [cosd(theta).*radius+centerpt(1)  sind(theta).*radius+centerpt(2)];
    oppositePts = [cosd(theta-180).*radius+centerpt(1)  sind(theta-180).*radius+centerpt(2)];
    if doplot
        plot(pts2scan(:,1),pts2scan(:,2),'-r',pts2scan(:,1),pts2scan(:,2),'ok')
    end
end
% --
if doplot
    plot( centerpt(1),centerpt(2),'or','MarkerFaceColor','r' )
%     plot( normalizedVirtualPolyCenter(1),normalizedVirtualPolyCenter(2),'ob','MarkerFaceColor','b' )
end
% --
allDist = nan(size(pts2scan,1),2);
for i=1:size(pts2scan,1)
%     plot( [pts2scan(i,1) oppositePts(i,1)] , [pts2scan(i,2) oppositePts(i,2)] ,':r' )
    X1 = [pts2scan(i,1) oppositePts(i,1)];   Y1 = [pts2scan(i,2) oppositePts(i,2)];
    iLineKeep = [];
    for iLine = 1:size(normalizedContourTemp,1)-1
        X2 = [normalizedContourTemp(iLine,1) normalizedContourTemp(iLine+1,1)];   Y2 = [normalizedContourTemp(iLine,2) normalizedContourTemp(iLine+1,2)];
        [x0,y0] = intersections(X1,Y1,X2,Y2);
        if ~isempty(x0) && ~isempty(y0)
            iLineKeep(end+1,:) = [ iLine x0 y0];
        end
    end
    if ~isempty(iLineKeep)
        % keep only shortest distance as the other intersection is an opposite line
        [minDist indMin] = min( sqrt( (pts2scan(i,1)-iLineKeep(:,2)).^2  +  (pts2scan(i,2)-iLineKeep(:,3)).^2 ) );
        coordIntersect = iLineKeep(indMin,2:3);
        if doplot
            plot( [pts2scan(i,1) coordIntersect(1)] , [pts2scan(i,2) coordIntersect(2)] ,'o-k' )
        end
        allDist(i,:) = [theta(i) minDist];
    end
end
% --
if doplot
    figure; plot(allDist(:,1),allDist(:,2),'o'); xlabel('\theta'); ylabel('normalized distance');
end

%%
function [HPZmask prop] = CalculateHPZmask( currentContour,HPZmask,dxmask,xmask,ymask,HPZareas )
in2 = false(size(HPZmask));
minXmask = max([min(currentContour(:,1))-dxmask*2  0]);
maxXmask = min([max(currentContour(:,1))+dxmask*2  max(xmask(:))]);
ind_j = nearestpoint( [minXmask maxXmask],xmask(1,:) );
minYmask = max([min(currentContour(:,2))-dymask*2  0]);
maxYmask = min([max(currentContour(:,2))+dymask*2  max(ymask(:))]);
ind_i = nearestpoint( [minYmask maxYmask],ymask(:,1) );
for i=ind_i(1):ind_i(2)
    for j=ind_j(1):ind_j(2)
        in2(i,j) = inpolygon(xmask(i,j),ymask(i,j),currentContour(:,1),currentContour(:,2));
    end
end
HPZmask(in2)                 = numel(HPZareas); % will set a unique number to that HPZ (i.e., from 1 to n)
% --
temp = zeros(size(HPZmask)); temp(in2) = 1;
prop = regionprops(temp,'area',...          % Scalar; the actual number of pixels in the region.
                        'BoundingBox',...   % The smallest rectangle containing the region
                        'Centroid',...      % 1-by-Q vector that specifies the center of mass of the region
                        'Eccentricity',...  % Scalar that specifies the eccentricity of the ellipse that 
                        ...                 % has the same second-moments as the region. The eccentricity 
                        ...                 % is the ratio of the distance between the foci of the ellipse and 
                        ...                 % its major axis length. The value is between 0 and 1. (0 and 1 are 
                        ...                 % degenerate cases; an ellipse whose eccentricity is 0 is actually a circle, 
                        ...                 % while an ellipse whose eccentricity is 1 is a line segment.)
                        'EquivDiameter',... % Scalar that specifies the diameter of a circle with the same area as the region
                        'EulerNumber',...   % Scalar that specifies the number of objects in the region minus the number of holes 
                        ...                 % in those objects.
                        'Extent',...        % Scalar that specifies the ratio of pixels in the region to pixels in the total bounding box
                        'MajorAxisLength',...% Scalar specifying the length (in pixels) of the major axis of the ellipse that has the same 
                        ...                  % normalized second central moments as the region.
                        'MinorAxisLength',...% Scalar; the length (in pixels) of the minor axis of the ellipse that has the same 
                        ...                  % normalized second central moments as the region.
                        'Orientation',...   % Scalar; the angle (in degrees ranging from -90 to 90 degrees) between the x-axis and the 
                        ...                 % major axis of the ellipse that has the same second-moments as the region.
                        'Perimeter',...     % Scalar; the distance around the boundary of the region. regionprops computes the perimeter 
                        ...                 % by calculating the distance between each adjoining pair of pixels around the border of the 
                        ...                 % region. If the image contains discontiguous regions, regionprops returns unexpected results.
                        'Solidity' );     	% Scalar specifying the proportion of the pixels in the convex hull that are also in the 
                                            % region. Computed as Area/ConvexArea.);
% -- 

%%
function irregularPolyCenter = my_irregularPolyCenter(x,y)
% NO GOOD AS FAR AS I'M CONCERNED
% --
fig1=figure('Position',[-1612 428 560 420]); hold on; plot(x,y,'o-');
% -- 
ok = 0;
xpoly = x;
ypoly = y;
while ~ok && (numel(x)>=2 && numel(y)>=2)
    min_max_x = [min(x) max(x)];
    min_max_y = [min(y) max(y)];
    center2try = [mean(min_max_x) mean(min_max_y)];
    hPlot = plot(center2try(1),center2try(2),'or');
    if ~inpolygon(center2try(1),center2try(2),xpoly,ypoly)
        x( x==min_max_x(1) ) = [];  x( x==min_max_x(2) ) = [];
        y( y==min_max_y(1) ) = [];  y( y==min_max_y(2) ) = [];
    else
        ok = 1;
    end
end
plot(center2try(1),center2try(2),'og','MarkerFaceColor','g');
irregularPolyCenter = center2try;


