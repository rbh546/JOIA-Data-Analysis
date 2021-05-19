function varargout = JOIAtool_2019(varargin)
%JOIATOOL is a GUI tool to analyze JOIA data.
%   JOIATOOL opens the GUI tool to analyze available JOIA data 
% 
%   Ridwan Hossain
%   Memorial University, St. John's, NL CANADA




%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JOIAtool_OpeningFcn, ...
                   'gui_OutputFcn',  @JOIAtool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
%     modified_gui_mainfcn(gui_State, varargin{:});
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%% --- Executes just before JOIAtool is made visible.
function JOIAtool_OpeningFcn(hObject, eventdata, handles, varargin)
% --
handles.output = hObject;
% -- INIT.
handles.options.grayShade2use = [0.7 0.7 0.7]; % [0 127 0]./255
handles.options.colorOrder = [0 180 0; 255 0 0; 191 191 0; 0 0 255]/255;
linkaxes([handles.ax2 handles.ax3 handles.ax4],'x')
% -- INIT. ...
widthTemp = 0.06; deltaY = -0.012;
    posTemp = get(handles.sensor1Check,'Position'); newPos = [posTemp(1)+posTemp(3) posTemp(2)+deltaY widthTemp posTemp(4)];
uicontrol(handles.uipanel12,'Style','text','Units','Normalized','Position',newPos,'String','(-)','ForegroundColor',handles.options.colorOrder(1,:),'BackgroundColor',get(handles.sensor1Check,'BackgroundColor'));
    posTemp = get(handles.sensor2Check,'Position'); newPos = [posTemp(1)+posTemp(3) posTemp(2)+deltaY widthTemp posTemp(4)];
uicontrol(handles.uipanel12,'Style','text','Units','Normalized','Position',newPos,'String','(-)','ForegroundColor',handles.options.colorOrder(2,:),'BackgroundColor',get(handles.sensor2Check,'BackgroundColor'));
    posTemp = get(handles.sensor3Check,'Position'); newPos = [posTemp(1)+posTemp(3) posTemp(2)+deltaY widthTemp posTemp(4)];
uicontrol(handles.uipanel12,'Style','text','Units','Normalized','Position',newPos,'String','(-)','ForegroundColor',handles.options.colorOrder(3,:),'BackgroundColor',get(handles.sensor3Check,'BackgroundColor'));
    posTemp = get(handles.sensor4Check,'Position'); newPos = [posTemp(1)+posTemp(3) posTemp(2)+deltaY widthTemp posTemp(4)];
uicontrol(handles.uipanel12,'Style','text','Units','Normalized','Position',newPos,'String','(-)','ForegroundColor',handles.options.colorOrder(4,:),'BackgroundColor',get(handles.sensor4Check,'BackgroundColor'));
% -- INIT. ...
handles.options.isFirstClick  = 1;
handles.options.isFirstLoad = 1;
% --
guidata(hObject, handles);
% -- Set objs. on/off until a file is loaded
SetObjOnOff(handles,'off','initial/load')
% --
set(handles.MainFig,'position',[-1670  40  3335  978])
% --

%% --- Outputs from this function are returned to the command line.
function varargout = JOIAtool_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%% --- Executes when user attempts to close MainFig.
function MainFig_CloseRequestFcn(hObject, eventdata, handles) %#ok<*INUSD>
delete(hObject);

%% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

%% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% --- Executes on selection change in ListFolders.
function ListFolders_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function ListFolders_CreateFcn(hObject, eventdata, handles)
handles = guihandles(hObject);
% --
handles.options.MainFolder = [cd,'\JOIA\Data\Data_1998'];
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --
guidata(hObject,handles);
% --

%% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles) %#ok<*INUSL>
MainFig_CloseRequestFcn(handles.MainFig,eventdata,guidata(handles.MainFig))
% --

%% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
ListFolders = get(handles.ListFolders,'String');
Val = get(handles.ListFolders,'Value');
currentFolder = ListFolders{Val}(1:8);
if ~strcmpi(currentFolder(1:2),'98')
    return;
end
handles.options.currentFolder = currentFolder;
handles.options.isFirstClick = 1;
% --

% load data
Data = readJOIAdata( [handles.options.MainFolder,filesep,currentFolder] );
    avgWindowSize = str2double(get(handles.avgWindowSize,'String'));
    if isnan(avgWindowSize) || isempty(avgWindowSize), avgWindowSize = 10; end
Data = performCalculations( Data , avgWindowSize );
set(handles.numFramesStatic,'String',[get(handles.numFramesStatic,'UserData'),num2str(numel(Data.calculations.timeVec))]);
% --
% set folder name
set(handles.currentFile,'String',['Current Data: ',currentFolder])
% --
% show the data
[Data, handles] = ShowAx1(Data,handles);
[Data, handles] = ShowAx2(Data,handles);
[Data, handles] = ShowAx3(Data,handles);
[Data, handles] = ShowAx4(Data,handles);
[Data, handles] = ShowAx5(Data,handles);
[Data, handles] = ShowAx6(Data,handles);
    Contour = ones(1,3)*1e6;  Peaks = ones(1,3)*2e6;
[Data, handles] = ShowAx10(Data,handles,Contour,Peaks);
[Data, handles] = ShowAx8(Data,handles);
[Data, handles, info] = ShowAx9(Data,handles);
% --
handles.options.isFirstLoad = 0;
handles.currentData = Data;
guidata(hObject,handles);
% -- show info on time series
set(handles.indentSpeed,    'String',[get(handles.indentSpeed,'UserData'),   ' ',info.IndentationSpeed ]);
set(handles.iceThickness,   'String',[get(handles.iceThickness,'UserData'),  ' ',info.IceThickness ]);
set(handles.iceStrength,    'String',[get(handles.iceStrength,'UserData'),   ' ',info.IceStrength ]);
set(handles.salinity,       'String',[get(handles.salinity,'UserData'),      ' ',info.Salinity ]);
set(handles.density,        'String',[get(handles.density,'UserData'),       ' ',info.Density ]);
set(handles.iceTemp,        'String',[get(handles.iceTemp,'UserData'),       ' ',info.IceTemperature]);
set(handles.iceType,        'String',[get(handles.iceType,'UserData'),       ' ',info.IceType ]);
set(handles.intactIce,      'String',[get(handles.intactIce,'UserData'),   	 ' ',info.IntactIceCondition ]);
set(handles.IceFailureMode, 'String',[get(handles.IceFailureMode,'UserData'),' ',info.IceFailureMode ]);
% --
iEdit = str2double(get(handles.iEdit,'String'));
if isempty(iEdit) || isnan(iEdit) || iEdit>numel(Data.calculations.timeVec), iEdit = 1; end
set(handles.iEdit,'String',num2str(iEdit));
if iEdit>1
    handles.options.isFirstClick    = 0;
    handles.options.isRestart       = 4;
    handles.options.isDoOneFrame    = true;
    guidata(hObject,handles);
    setappdata(handles.MainFig,'isStoppedFlag',0)
    GetRolling(hObject, eventdata, handles)
    handles.options.isFirstClick    = 1;
    guidata(hObject,handles);
end
% -- Set objs. on/off when a file is loaded
SetObjOnOff(handles,'on','initial/load')
% -- End loading data
 
%% Perform Calculations
function Data = performCalculations( Data,avgWindowSize )
% -- calculate forces, pressures and areas, based on tactile data to ensure matches the ASG file provided
[TactileCalculatedForce, area, press, ForcesPerSensor, AreasPerSensor, PressPerSensor] = doCalculations(Data);
% -- calculate global pressures from data in ASG files
Press_ASG = (Data.TactileSensor.ASGfile.force.Force./1000) ./ (Data.TactileSensor.ASGfile.area.Area*1e-6); % converted from mm^2 to m^2 to get kPa, then /1000 to get MPa
% -- apply a smoothing avg. filter on the global pressures and areas calculated based on ASF file (movie file)
press_smooth = smoothData( press,avgWindowSize );
area_smooth  = smoothData( area,avgWindowSize );
% -- calculate back a smoothed version of force
TactileCalculatedForce_smooth = press_smooth.*area_smooth*1000;
% -- calculate the slope of pressure curve
diff_press          = diff(press);        	diff_press          = [nan;diff_press];
diff_press_smooth   = diff(press_smooth);   diff_press_smooth   = [nan;diff_press_smooth];
% -- calculate the slope of force curve
diff_force          = diff(TactileCalculatedForce);        	diff_force          = [nan;diff_force];
diff_force_smooth   = diff(TactileCalculatedForce_smooth); 	diff_force_smooth   = [nan;diff_force_smooth];
% --
maxV = nan(size(Data.TactileSensor.ASFfile.Pressures,3),1);
for i=1:size(Data.TactileSensor.ASFfile.Pressures,3), maxV(i) = max(max(Data.TactileSensor.ASFfile.Pressures(:,:,i))); end
maxV = max(maxV);
satPressFlag = 0;
if abs(maxV-Data.TactileSensor.ASFfile.SATURATION_PRESSURE)<0.01, satPressFlag = 1; end
% --
TotalAreaAllPanels = Data.TactileSensor.ASFfile.ROWS * Data.TactileSensor.ASFfile.ROW_SPACING * Data.TactileSensor.ASFfile.COLS * Data.TactileSensor.ASFfile.COL_SPACING *1e-6; % convert from mm^2 to m^2;
% -- smooth data by sensor
PressPerSensor_smooth = smoothData( PressPerSensor,avgWindowSize );
ForcesPerSensor_smooth = smoothData( ForcesPerSensor,avgWindowSize );
AreasPerSensor_smooth = smoothData( AreasPerSensor,avgWindowSize );
% --
% look how many cells and fit to screen
[m n] = size(Data.TactileSensor.ASFfile.Pressures(:,:,1));
if m~=44, error('Have a look, # of rows is different than 44'); return; end
if n~=176
    filledPortion = (176-n)/2;
    panelPos = [filledPortion+1  filledPortion+n];
    dataTemp = zeros(44,176,size(Data.TactileSensor.ASFfile.Pressures,3));
    for i=1:size(Data.TactileSensor.ASFfile.Pressures,3)
        dataTemp(:,panelPos(1):panelPos(2),i) = Data.TactileSensor.ASFfile.Pressures(:,:,i);
    end
    Data.TactileSensor.ASFfile.Pressures = dataTemp;
    Data.TactileSensor.ASFfile.COLS = 176;
end
% -- Create X and Y matrices + time vector
x = 0:Data.TactileSensor.ASFfile.COL_SPACING:Data.TactileSensor.ASFfile.COL_SPACING*Data.TactileSensor.ASFfile.COLS-Data.TactileSensor.ASFfile.COL_SPACING;
y = 0:Data.TactileSensor.ASFfile.ROW_SPACING:Data.TactileSensor.ASFfile.ROW_SPACING*Data.TactileSensor.ASFfile.ROWS-Data.TactileSensor.ASFfile.ROW_SPACING;
[X Y] = meshgrid(x,y);
timeVec = 0:Data.TactileSensor.ASFfile.SECONDS_PER_FRAME:Data.TactileSensor.ASFfile.SECONDS_PER_FRAME*(Data.TactileSensor.ASFfile.END_FRAME-Data.TactileSensor.ASFfile.START_FRAME);
% -- options for plotting
fraction = 0.2;
minX = -fraction*Data.TactileSensor.ASFfile.COL_SPACING;
maxX = max(x)+fraction*Data.TactileSensor.ASFfile.COL_SPACING;
minY = -fraction*Data.TactileSensor.ASFfile.ROW_SPACING;
maxY = max(y)+fraction*Data.TactileSensor.ASFfile.ROW_SPACING;

% -- KEEP all in Data
calculations.TactileCalculatedForce = TactileCalculatedForce;
calculations.area = area;
calculations.press = press;
calculations.ForcesPerSensor = ForcesPerSensor;
calculations.AreasPerSensor = AreasPerSensor;
calculations.PressPerSensor = PressPerSensor;
calculations.Press_ASG = Press_ASG;
calculations.press_smooth = press_smooth;
calculations.area_smooth = area_smooth;
calculations.TactileCalculatedForce_smooth = TactileCalculatedForce_smooth;
calculations.PressPerSensor_smooth = PressPerSensor_smooth;
calculations.ForcesPerSensor_smooth = ForcesPerSensor_smooth;
calculations.AreasPerSensor_smooth = AreasPerSensor_smooth;
calculations.diff_press = diff_press;
calculations.diff_press_smooth = diff_press_smooth;
calculations.diff_force = diff_force;
calculations.diff_force_smooth = diff_force_smooth;
calculations.TotalAreaAllPanels = TotalAreaAllPanels;
calculations.satPressFlag = satPressFlag;
calculations.x = x;
calculations.y = y;
calculations.X = X;
calculations.Y = Y;
calculations.timeVec = timeVec;
calculations.minX = minX;
calculations.maxX = maxX;
calculations.minY = minY;
calculations.maxY = maxY;
% --
Data.calculations = calculations;
% -- End of calculation

%% ShowAx1 - Plots Tactile pressure sensor data on Ax1
function [Data, handles] = ShowAx1(Data,handles)
% --
minX = Data.calculations.minX;
minY = Data.calculations.minY;
maxX = Data.calculations.maxX;
maxY = Data.calculations.maxY;
X = Data.calculations.X;
Y = Data.calculations.Y;
x = Data.calculations.x;
y = Data.calculations.y;
satPressFlag = Data.calculations.satPressFlag;
% --
% look how many cells and fit to screen
[m n] = size(Data.TactileSensor.ASFfile.Pressures(:,:,1));
if m~=44, error('Have a look, # of rows is different than 44'); return; end
if n~=176
    filledPortion = (176-n)/2;
    panelPos = [filledPortion+1  filledPortion+n];
    dataTemp = zeros(44,176,size(Data.TactileSensor.ASFfile.Pressures,3));
    for i=1:size(Data.TactileSensor.ASFfile.Pressures,3)
        dataTemp(:,panelPos(1):panelPos(2),i) = Data.TactileSensor.ASFfile.Pressures(:,:,i);
    end
    Data.TactileSensor.ASFfile.Pressures = dataTemp;
    Data.TactileSensor.ASFfile.COLS = 176;
else
    panelPos = [ 45 89 133];
end
% --
set(handles.ax1,...
    'XLim',[minX-3*Data.TactileSensor.ASFfile.COL_SPACING   maxX+1*Data.TactileSensor.ASFfile.COL_SPACING],...
    'YLim',[minY-3*Data.TactileSensor.ASFfile.ROW_SPACING   maxY+1*Data.TactileSensor.ASFfile.ROW_SPACING] )
axes( handles.ax1 ); %#ok<*MAXES>
hold(handles.ax1,'on'); axis equal
% --

% --
if ~handles.options.isFirstLoad
    % if not first load of the session, do not recreate plot objects on top 
    % of existing ones - erase old ones first
    delete(handles.currentData.plotHandles.hSurf)
    delete(handles.currentData.plotHandles.hCbar1)
%     delete(handles.currentData.plotHandles.hCbarTitle1)
    try delete(handles.currentData.plotHandles.satPress); end
end
% --

% --
hSurf = surf( X,Y,Data.TactileSensor.ASFfile.Pressures(:,:,1) );
set(hSurf,'Tag','hSurf','CDataMode','auto');
view(2); shading interp;
% --
ColorMap = GetColorMap(1);
colormap(handles.ax1,ColorMap)
colorbar;
caxis( handles.ax1, [ 0  max(Data.TactileSensor.ASFfile.Pressures(:)) ] )  % [ min(Data.TactileSensor.ASFfile.Pressures(:))  max(Data.TactileSensor.ASFfile.Pressures(:)) ]
% hCbar1=handles.ax1;
hCbar1 = cbfreeze(handles.ax1); % freeze colorbar
% hCbarTitle1 = title(hCbar1,['Local Pressure (',Data.TactileSensor.ASFfile.UNITS,')']);
hCbarTitle1 = cblabel(['Local Pressure (',Data.TactileSensor.ASFfile.UNITS,')'],'Rotation',0,'units','normalized','Position',[0.5 1.07 0]);
% --
for j=1:10:numel(x), plot( [x(j) x(j)],[min(y) max(y)] ,':','Color',[0.8 0.8 0.8] ); end
for j=1:5:numel(y), plot( [min(x) max(x)],[y(j) y(j)] ,':','Color',[0.8 0.8 0.8] ); end
plot( [minX maxX maxX minX minX],[minY minY maxY maxY minY] ,'-k' );
% --
panelPosHdl(1,1) = plot( [x(panelPos(1)-1) x(panelPos(1)-1)],[min(y) max(y)] ,'--','Color',[0.6 0.6 0.6],'LineWidth',1,'Visible','off','Tag','panelPosHdl' );
panelPosHdl(end+1,1) = plot( [x(panelPos(2)+1) x(panelPos(2)+1)],[min(y) max(y)] ,'--','Color',[0.6 0.6 0.6],'LineWidth',1,'Visible','off','Tag','panelPosHdl' );
if numel(panelPos)<3
    set(panelPosHdl,'Visible','on');
else
    panelPosHdl(end+1,1) = plot( [x(panelPos(3)+1) x(panelPos(3)+1)],[min(y) max(y)] ,'--','Color',[0.6 0.6 0.6],'LineWidth',1,'Visible','off','Tag','panelPosHdl' );
end
% --
xlabel( ['Horizontal Distance (',Data.TactileSensor.ASFfile.COL_SPACING_units,')'] )
ylabel( ['Vertical Distance (',Data.TactileSensor.ASFfile.ROW_SPACING_units,')'] )
% --
satPress = [];
hold(handles.ax1,'off');
% -- freeze colormap
freezeColors(handles.ax1); % freeze colormap
% -- KEEP all in Data
plotHandles.panelPos = panelPos;
plotHandles.hSurf = hSurf;
%plotHandles.hCbar1 = hCbar1;
plotHandles.hCbarTitle1 = hCbarTitle1;
plotHandles.satPress = satPress;
plotHandles.panelPosHdl = panelPosHdl;
% --
Data.plotHandles = plotHandles;
% -- End of function for Ax1
    
%% ShowAx2 - Plots Forces on Ax2
function [Data, handles] = ShowAx2(Data,handles)
% --
plotSmoothFlag = get(handles.useSmoothData,'Value');
timeVec = Data.calculations.timeVec;
ForcesPerSensor = Data.calculations.ForcesPerSensor;
ForcesPerSensor_smooth = Data.calculations.ForcesPerSensor_smooth;
TactileCalculatedForce = Data.calculations.TactileCalculatedForce;
TactileCalculatedForce_smooth = Data.calculations.TactileCalculatedForce_smooth;
grayShade2use = handles.options.grayShade2use;
colorOrder = handles.options.colorOrder;
% --
axes( handles.ax2 ); %#ok<*MAXES>
hold(handles.ax2,'on'); grid on;
set(handles.ax2,'YLimMode','auto');
% --

% --
if ~handles.options.isFirstLoad
    % if not first load of the session, do not recreate plot objects on top 
    % of existing ones - erase old ones first
    delete(handles.currentData.plotHandles.hLjoia1)
    delete(handles.currentData.plotHandles.hLjoia12)
    delete(handles.currentData.plotHandles.hLjoia13)
    delete(handles.currentData.plotHandles.hLsensor(:,1))
    delete(handles.currentData.plotHandles.hL2)
end
% --

% --
hLjoia1 = plot( Data.TactileSensor.ASGfile.force.Time , Data.TactileSensor.ASGfile.force.Force ,'--','LineWidth',1,'Color',[.7 .7 .7] );
set(hLjoia1,'Tag','hLjoia1','Visible','off');
% --
if plotSmoothFlag, finalData = ForcesPerSensor_smooth; else finalData = ForcesPerSensor; end %#ok<*UNRCH>
for iSens = 1:4
    hLsensor(iSens,1) = plot( timeVec , finalData(:,iSens) ,'-','Color',colorOrder(iSens,:),'LineWidth',1 ,'Tag','hLsensor1','Visible','off'); %#ok<*AGROW>
end
% --
ylabel( ['Total Force (kN)'] ) %#ok<*NBRAK>
set( handles.ax2,'XLim',[0 ceil(max(timeVec))] ,'XTickLabel','' )
% -- show calculated forces based on tactile data to ensure matches the ASG file provided
if plotSmoothFlag, finalData = TactileCalculatedForce_smooth; else finalData = TactileCalculatedForce; end %#ok<*UNRCH>
% -- show
hLjoia13 = []; if plotSmoothFlag,  hLjoia13 = plot( timeVec , TactileCalculatedForce ,'-','Color',grayShade2use-0.2,'LineWidth',1 );  end
hLjoia12 = plot( timeVec , finalData ,'-k','LineWidth',2 );
% --
ax2YLim = get(handles.ax2,'YLim'); set(handles.ax2,'YLim',ax2YLim);
hL2 = plot( [timeVec(1) timeVec(1)],ax2YLim , 'r-','LineWidth',2 );
set(hL2,'Tag','hL2');
% --
hold(handles.ax2,'off');
% -- KEEP all in Data
plotHandles = Data.plotHandles;
% --
plotHandles.hLjoia1 = hLjoia1;
plotHandles.hLjoia12 = hLjoia12;
plotHandles.hLjoia13 = hLjoia13;
plotHandles.hLsensor = hLsensor;
plotHandles.hL2 = hL2;
plotHandles.origPosAx2 = get(handles.ax2,'Position');
% --
Data.plotHandles = plotHandles;
% -- End of function for Ax2

%% ShowAx3 - Plot Areas on Ax3
function [Data, handles] = ShowAx3(Data,handles)
% --
plotSmoothFlag = get(handles.useSmoothData,'Value');
timeVec = Data.calculations.timeVec;
AreasPerSensor = Data.calculations.AreasPerSensor;
AreasPerSensor_smooth = Data.calculations.AreasPerSensor_smooth;
area = Data.calculations.area;
area_smooth = Data.calculations.area_smooth;
TotalAreaAllPanels = Data.calculations.TotalAreaAllPanels;
grayShade2use = handles.options.grayShade2use;
colorOrder = handles.options.colorOrder;
hLsensor = Data.plotHandles.hLsensor;
% --
axes( handles.ax3 ); %#ok<*MAXES>
hold(handles.ax3,'on'); grid on;
set(handles.ax3,'YLimMode','auto');
% --

% --
if ~handles.options.isFirstLoad
    % if not first load of the session, do not recreate plot objects on top 
    % of existing ones - erase old ones first
    delete(handles.currentData.plotHandles.hLjoia2)
    delete(handles.currentData.plotHandles.hLjoia22)
    delete(handles.currentData.plotHandles.hLjoia23)
    delete(handles.currentData.plotHandles.hLsensor(:,2))
    delete(handles.currentData.plotHandles.hL3)
    delete(handles.currentData.plotHandles.ax3R)
end
% --

% --
hLjoia2 = plot( Data.TactileSensor.ASGfile.area.Time , (Data.TactileSensor.ASGfile.area.Area*1e-6) ,'--','LineWidth',1,'Color',[.7 .7 .7] );
set(hLjoia2,'Tag','hLjoia2','Visible','off');
% --
if plotSmoothFlag, finalData = AreasPerSensor_smooth; else finalData = AreasPerSensor; end %#ok<*UNRCH>
for iSens = 1:4
    hLsensor(iSens,2) = plot( timeVec , finalData(:,iSens) ,'-','Color',colorOrder(iSens,:),'LineWidth',1 ,'Tag','hLsensor2','Visible','off');
end
% --
ylabel( ['Contact Area (m^2)'] )
set( handles.ax3,'XLim',[0 ceil(max(timeVec))] ,'XTickLabel','' )
% -- show calculated forces based on tactile data to ensure matches the ASG file provided
if plotSmoothFlag, finalDataArea = area_smooth; else finalDataArea = area; end %#ok<*UNRCH>
hLjoia23 = []; if plotSmoothFlag,  hLjoia23 = plot( timeVec , area ,'-','Color',grayShade2use-0.2,'LineWidth',1 );  end
hLjoia22 = plot( timeVec , finalDataArea ,'-k','LineWidth',2 );
% --
ax3YLim = get(handles.ax3,'YLim'); set(handles.ax3,'YLim',ax3YLim);
hL3 = plot( [timeVec(1) timeVec(1)],ax3YLim , 'r-','LineWidth',2 );
set(hL3,'Tag','hL3');

hold(handles.ax3,'off');
% -- KEEP all in Data
plotHandles = Data.plotHandles;
% --
plotHandles.hLjoia2 = hLjoia2;
plotHandles.hLjoia22 = hLjoia22;
plotHandles.hLjoia23 = hLjoia23;
plotHandles.hLsensor = hLsensor;
% plotHandles.ax3R = ax3R;
plotHandles.hL3 = hL3;
plotHandles.origPosAx3 = get(handles.ax3,'Position');
% --
Data.plotHandles = plotHandles;
% -- End of function for Ax3

%% ShowAx4 - Plots Pressures on Ax4
function [Data, handles] = ShowAx4(Data,handles)
% --
plotSmoothFlag = get(handles.useSmoothData,'Value');
timeVec = Data.calculations.timeVec;
Press_ASG = Data.calculations.Press_ASG;
PressPerSensor = Data.calculations.PressPerSensor;
PressPerSensor_smooth = Data.calculations.PressPerSensor_smooth;
press_smooth = Data.calculations.press_smooth;
press = Data.calculations.press;
grayShade2use = handles.options.grayShade2use;
colorOrder = handles.options.colorOrder;
hLsensor = Data.plotHandles.hLsensor;
% --
axes( handles.ax4 ); %#ok<*MAXES>
hold(handles.ax4,'on'); grid on;
set(handles.ax4,'YLimMode','auto');
% --

% --
if ~handles.options.isFirstLoad
    % if not first load of the session, do not recreate plot objects on top 
    % of existing ones - erase old ones first
    delete(handles.currentData.plotHandles.hLjoia3)
    delete(handles.currentData.plotHandles.hLjoia32)
    delete(handles.currentData.plotHandles.hLjoia33)
    delete(handles.currentData.plotHandles.hLsensor(:,3))
    delete(handles.currentData.plotHandles.hL4)
end
% --

% --
hLjoia3 = plot( timeVec , Press_ASG ,'--','LineWidth',1,'Color',[.7 .7 .7] );
set(hLjoia3,'Tag','hLjoia3','Visible','off');
% --
if plotSmoothFlag, finalData = PressPerSensor_smooth; else finalData = PressPerSensor; end %#ok<*UNRCH>
for iSens = 1:4
    hLsensor(iSens,3) = plot( timeVec , finalData(:,iSens) ,'-','Color',colorOrder(iSens,:),'LineWidth',1 ,'Tag','hLsensor3','Visible','off');
end
% --
ylabel( ['Mean Contact Pressure (MPa)'] )
set( handles.ax4,'XLim',[0 ceil(max(timeVec))] )
% -- calculate global pressures based on tactile data to ensure matches the ASG file provided
if plotSmoothFlag, finalData = press_smooth; else finalData = press; end %#ok<*UNRCH>
% -- show
hLjoia33 = []; if plotSmoothFlag,  hLjoia33 = plot( timeVec , press ,'-','Color',grayShade2use-0.2,'LineWidth',1 );  end
hLjoia32 = plot( timeVec , finalData ,'-k','LineWidth',2 );
% --
ax4YLim = get(handles.ax4,'YLim'); set(handles.ax4,'YLim',ax4YLim);
hL4 = plot( [timeVec(1) timeVec(1)],ax4YLim , 'r-','LineWidth',2 );
set(hL4,'Tag','hL4');
% --
hold(handles.ax4,'off');
% -- KEEP all in Data
plotHandles = Data.plotHandles;
% --
plotHandles.hLjoia3 = hLjoia3;
plotHandles.hLjoia32 = hLjoia32;
plotHandles.hLjoia33 = hLjoia33;
plotHandles.hLsensor = hLsensor;
plotHandles.hL4 = hL4;
plotHandles.origPosAx4 = get(handles.ax4,'Position');
% --
Data.plotHandles = plotHandles;
% -- End of function for Ax4

%% ShowAx5 - Plots pressure-area cloud on Ax5
function [Data, handles] = ShowAx5(Data,handles)
% --
plotSmoothFlag = get(handles.useSmoothData,'Value');
area = Data.calculations.area;
press = Data.calculations.press;
area_smooth = Data.calculations.area_smooth;
press_smooth = Data.calculations.press_smooth;
grayShade2use = handles.options.grayShade2use;
% --
axes( handles.ax5 ); %#ok<*MAXES>
hold(handles.ax5,'on'); grid on;
set(handles.ax5,'YLimMode','auto');
% --

% --
if ~handles.options.isFirstLoad
    % if not first load of the session, do not recreate plot objects on top 
    % of existing ones - erase old ones first
    delete(handles.currentData.plotHandles.hL5)
    delete(handles.currentData.plotHandles.hL52)
    delete(handles.currentData.plotHandles.hL53)
    delete(handles.currentData.plotHandles.hL5zoom)
    delete(handles.currentData.plotHandles.hL5_trace)
end
% --

% --
if plotSmoothFlag, finalData_press = [area_smooth,press_smooth]; else finalData_press = [area,press]; end
hL53 = [];
if plotSmoothFlag
    hL53 = plot( area,press ,'o','MarkerEdgeColor',grayShade2use,'MarkerFaceColor',grayShade2use,'MarkerSize',2 );
end
hL52 = plot( finalData_press(:,1),finalData_press(:,2) ,'o','MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',3 );
hL5zoom = plot( nan,nan ,'ok','MarkerEdgeColor',[0.101960784313725  0.250980392156863  0.4],'MarkerSize',4 ); set(hL5zoom,'Tag','hL5zoom','UserData',{finalData_press(:,1),finalData_press(:,2)});
% --
xlabel( ['Contact Area (m^2)'] )
ylabel( ['Mean Contact Pressure (MPa)'] )
XLim = get(handles.ax3,'YLim');  YLim = get(handles.ax4,'YLim');
set(handles.ax5,'XLim',XLim,'YLim',YLim);
% --
hL5 = plot( finalData_press(1,1),finalData_press(1,2) ,'or','MarkerFaceColor','r','MarkerSize',4 );
set(hL5,'Tag','hL5');
hL5_trace = plot( nan,nan,'-','LineWidth',1,'Color',[177 14 14]/255);
set(hL5_trace,'Tag','hL5_trace');
% --
hold(handles.ax5,'off');
% -- KEEP all in Data
plotHandles = Data.plotHandles;
% --
plotHandles.hL5 = hL5;
plotHandles.hL52 = hL52;
plotHandles.hL53 = hL53;
plotHandles.hL5zoom = hL5zoom;
plotHandles.hL5_trace = hL5_trace;
% --
Data.plotHandles = plotHandles;
% --


function update_display(cursorPositions,data,TotalAreaAllPanels)
xData = data{1};
% -- find graphic objects
ax2     = findobj('Tag','ax2');
ax3     = findobj('Tag','ax3');
ax3R    = findobj('Tag','ax3R');
ax4     = findobj('Tag','ax4');
EdStart = findobj('Tag','EdStart');
EdEnd   = findobj('Tag','EdEnd');
hL5zoom = findobj('Tag','hL5zoom');
hL6zoom = findobj('Tag','hL6zoom');
% -- calculate limits and corresponding indices, ...
I1 = find(cursorPositions(1) >= xData,1,'last');
I2 = find(cursorPositions(2) <= xData,1,'first');
if isempty(I1) || isnan(I1) || I1<1,            I1 = 1;            end;
if isempty(I2) || isnan(I2) || I2>numel(xData), I2 = numel(xData); end
cursorPositions(1) = xData(I1);
cursorPositions(2) = xData(I2);
% -- set new axes limits - 1: auto to get new limits, ...
set(ax2,'XLim',cursorPositions,'YLimMode','auto'); % note, ax3 and ax4 will follow the xlimas they are already linked together
set(ax3,'YLimMode','auto');
set(ax4,'YLimMode','auto');
% ... 2: fix them so that when running it does not account for parts that we don't see
YLim = get(ax2,'YLim'); set(ax2,'YLim',YLim); 
YLim = get(ax3,'YLim'); set(ax3,'YLim',YLim); 
YLim = get(ax4,'YLim'); set(ax4,'YLim',YLim); 
% -- Update times in the edit boxes
set(EdStart,'string',['min. Time = ',num2str(cursorPositions(1),'%0.2f'),' sec.  ;  (Frame #: ',num2str(I1),')'],'UserData',I1);
set(EdEnd  ,'string',['max. Time = ',num2str(cursorPositions(2),'%0.2f'),' sec.  ;  (Frame #: ',num2str(I2),')'],'UserData',I2);
% -- update right y-axis on area plot
ax3YTicks = get(ax3,'YTick');
ax3RYTicks = ax3YTicks/TotalAreaAllPanels;
tempLab(1,:) = '0   '; for iLab=2:numel(ax3YTicks), tempLab(iLab,:) = num2str( ax3RYTicks(iLab) ,'%0.2f' ); end
set(ax3R,'YTick',ax3RYTicks,'YLim',[0 max(ax3RYTicks)],'YTickLabel',tempLab);
% --
hL5UserData = get(hL5zoom,'UserData');
hL6UserData = get(hL6zoom,'UserData');
set(hL5zoom,'XData',hL5UserData{1}(I1:I2),'YData',hL5UserData{2}(I1:I2))
set(hL6zoom,'XData',hL6UserData{1}(I1:I2),'YData',hL6UserData{2}(I1:I2))
% --


%% GetRolling - Start showing in loop
function GetRolling(hObject, eventdata, handles)
% -- turn objects enabling off while rolling
SetObjOnOff(handles,'off')
% -- % retrieve important data
Data            = handles.currentData;
isRestart       = handles.options.isRestart;
isDoOneFrame    = handles.options.isDoOneFrame;
timeVec         = Data.calculations.timeVec;
currentFolder   = handles.options.currentFolder;
plotSmoothFlag  = get(handles.useSmoothData,'Value');
DoAx10Flag      = get(handles.EnableAx10,'Value');
if plotSmoothFlag
    finalData_press         = [Data.calculations.area_smooth , Data.calculations.press_smooth];
    finalData_diff_press    = [Data.calculations.area_smooth , Data.calculations.diff_press_smooth]; 
else
    finalData_press         = [Data.calculations.area , Data.calculations.press];
    finalData_diff_press    = [Data.calculations.area , Data.calculations.diff_press];
end
% --
ColorMap1 = GetColorMap(1);
ColorMap2 = GetColorMap(2);
% --
i = str2double(get(handles.iEdit,'String'));
if isRestart>0
    if      isRestart==2, i = i-1;
    elseif  isRestart==3, i = i+1;
    elseif  isRestart==4, % i=i ; i.e., stay the same
    end
end
% --
if (isRestart==0&&~isnumeric(i)) || isnan(i) || isempty(i) || i>numel(timeVec) || i==1
    i = 2;
end
% --
try delete( findobj(handles.MainFig,'Tag','hText') ); end
if ~DoAx10Flag
    % erase existing stuff on ax10 if not used
    try delete( findobj(handles.MainFig,'Tag','hSurf10') ); end
    try delete( findobj(handles.MainFig,'Tag','hPlot3_ax10') ); end
    try delete( findobj(handles.MainFig,'Tag','hdl_contour_ax10') ); end
        xlim = get(handles.ax10,'XLim');  ylim = get(handles.ax10,'YLim');
    text( mean(xlim),mean(ylim),'\bfNOTE:\rm OPTION TO SHOW HPZs CONTOURS AND PEAKS NOT ENABLED' ,...
          'Parent',handles.ax10,'HorizontalAlignment','center','FontWeight','bold','Tag','hText' );
end
% --
while 1  
    % --
    set(handles.frameStatic,'String',[get(handles.frameStatic,'UserData'),num2str(i,'%0.0f')])
    setappdata(handles.MainFig,'current_i',i)
    iPercentage = (i-1)/numel(timeVec)*100;
    % --
    if ~getappdata(handles.MainFig,'isStoppedFlag') % used to stop the loop
        title(handles.ax1, ['Time since beginning: ',num2str(timeVec(i-1),'%0.2f'),' sec. (',num2str(iPercentage,'%1.0f'),'% of data)   -   Folder: ',currentFolder] )
        % --
        currentData = Data.TactileSensor.ASFfile.Pressures(:,:,i);
        % --  show contours of HPZs
        if DoAx10Flag
            axes(handles.ax10) %#ok<LAXES>
            foundHPZs_data = findHPZs( Data.calculations.X,Data.calculations.Y,currentData,handles );
            peaks       = foundHPZs_data.peaks;
            c_contour	= foundHPZs_data.c_contour;
            ZLIM = get( handles.ax10 , 'ZLim' );
            peaks(:,3) = ZLIM(2)*20;
            c_contour(:,3) = ZLIM(2)*10;
            [Data, handles] = ShowAx10(Data,handles,c_contour,peaks);
        end
        % --  
        axes(handles.ax8) %#ok<LAXES>
        set( Data.plotHandles.hSurf,'CData',currentData,'ZData',currentData )
        % --  show dP/dt data and set colormaps
        if i>1
            previousData = Data.TactileSensor.ASFfile.Pressures(:,:,i-1);
            temp = currentData - previousData;
            set( Data.plotHandles.hSurfCmp,'CData',temp,'ZData',temp )
            colormap(handles.ax8,ColorMap2) % need to re-specify the colormap to avoid messing with other colormap in ax1
            freezeColors(handles.ax8); % need to freeze colormap because we will set a different one to ax1 below
        end
        % -- need to specify the colormap for ax1 and all other that aren't frozen (i.e., ax10) here as it has been lost above when mocking with ax8              
        colormap(handles.ax1,ColorMap1)          
        % -- vertical red lines
        set(Data.plotHandles.hL2,'XData',[timeVec(i) timeVec(i)]);
        set(Data.plotHandles.hL3,'XData',[timeVec(i) timeVec(i)]);
        set(Data.plotHandles.hL4,'XData',[timeVec(i) timeVec(i)]);
        % -- dots and their time traces
        iMin = max([1 i-25]);
        set(Data.plotHandles.hL5,       'XData',finalData_press(i,1),          'YData',finalData_press(i,2));
        set(Data.plotHandles.hL5_trace, 'XData',finalData_press(iMin:i,1),     'YData',finalData_press(iMin:i,2));
        set(Data.plotHandles.hL6,       'XData',finalData_diff_press(i,1),     'YData',finalData_diff_press(i,2));
        set(Data.plotHandles.hL6_trace, 'XData',finalData_diff_press(iMin:i,1),'YData',finalData_diff_press(iMin:i,2));
        % --
        drawnow
    else       
        break
    end
    % --
    i = i + 1;
    if i>numel(timeVec), i=2; end
    % --
    userSpeed = str2double(get(handles.speedEdit,'String'));
    if ~isempty(userSpeed) && isnumeric(userSpeed) && userSpeed>0 && userSpeed<100
        pause(1/userSpeed/5);
    end
    % --
    if isDoOneFrame, break; end
    % --
end
% --
i = getappdata(handles.MainFig,'current_i');
set(handles.iEdit,'String',num2str(i,'%0.0f'),'Enable','on')
% --
handles.options.isRestart = isRestart;
handles.options.isDoOneFrame = isDoOneFrame;
guidata(hObject,handles);
% -- turn objects enabling on
SetObjOnOff(handles,'on')
% -- End of GetRolling


%% A custom colormap was used for showing the tactile data on Ax1
function ColorMap = GetColorMap(ind)
% --
if ind==1
    inputColorMap = [ 1   255 255 255; % [IND R G B]
                      11  0   0   143;
                      18  0   0   255;
                      26  0   255 255;
                      35  255 255 0  ;
                      51  255 0   0  ;
                      86  65  0   0  ;
                     ];
    % adjust a little if wanted, comment newxt line if not
    inputColorMap(2:5,1) = inputColorMap(2:5,1) - 5;
elseif ind==2
    inputColorMap = cbrewer('div','RdBu',11);
    inputColorMap(6,:) = [1 1 1];
    inputColorMap = [[1:size(inputColorMap,1)]' inputColorMap*255];
elseif ind==3
    inputColorMap = [ 1   255 255 255; % [IND R G B]
                      25  0   0   143;
                      30  0   0   255;
                      35  0   255 255;
                      45  255 255 0  ;
                      65  255 0   0  ;
                      86  65  0   0  ;
                     ];
end
% --
if ind==3
    fullIND = linspace(0,1,max(inputColorMap(:,1)));
else
    fullIND = linspace(0,1,125);
end
tempIND = (inputColorMap(:,1)-1)/(inputColorMap(end,1)-1);
ColorMap = [interp1(tempIND,inputColorMap(:,2),fullIND)'   interp1(tempIND,inputColorMap(:,3),fullIND)'   interp1(tempIND,inputColorMap(:,4),fullIND)'];
ColorMap = ColorMap/max(ColorMap(:));
% -- End of colormap

%% --- Executes during object creation, after setting all properties.
function speedStatic_CreateFcn(hObject, eventdata, handles)

%%
function iEdit_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function iEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on button press in stopButton.
function stopButton_Callback(hObject, eventdata, handles)
% --
if handles.options.isFirstClick
    handles.options.isFirstClick    = 0;
    handles.options.isRestart       = 0;
    handles.options.isDoOneFrame    = false;
    setappdata(handles.MainFig,'isStoppedFlag',0)
    setappdata(handles.MainFig,'current_i',0)
else
    handles.options.isRestart    = 1;
    handles.options.isDoOneFrame = false;
    current = getappdata(handles.MainFig,'isStoppedFlag');
    setappdata(handles.MainFig,'isStoppedFlag',abs(current-1) )
end
guidata(hObject,handles);
% --
if ~getappdata(handles.MainFig,'isStoppedFlag')
	GetRolling(hObject, eventdata, handles)
end
% -- End of stopButton callback function

%% --- Executes on button press in goDnButton.
function goDnButton_Callback(hObject, eventdata, handles)
% --
setappdata(handles.MainFig,'isStoppedFlag',0 )
handles.options.isRestart 	 = 2;
handles.options.isDoOneFrame = true;
guidata(hObject,handles);
% --
GetRolling(hObject, eventdata, handles)
% -- End of goDnButton

%% --- Executes on button press in goUpButton.
function goUpButton_Callback(hObject, eventdata, handles)
% --
setappdata(handles.MainFig,'isStoppedFlag',0 )
handles.options.isRestart 	 = 3;
handles.options.isDoOneFrame = true;
guidata(hObject,handles);
% --
GetRolling(hObject, eventdata, handles)
% -- End of goUpButton

%% --- Executes on button press in joiaTimeTraceCheck.
function joiaTimeTraceCheck_Callback(hObject, eventdata, handles)
% --
vis = 'off'; 
if get(handles.joiaTimeTraceCheck,'Value'), vis = 'on'; end;
set(handles.currentData.plotHandles.hLjoia1,'Visible',vis);
if ~(strcmpi(vis,'on') && get(handles.forceOnlyCheck,'Value'))
    set(handles.currentData.plotHandles.hLjoia2,'Visible',vis);
    set(handles.currentData.plotHandles.hLjoia3,'Visible',vis);
end
% --
guidata(hObject,handles);
% --

%% --- Executes on button press in sensor1Check.
function sensor1Check_Callback(hObject, eventdata, handles)
setVisSensors(hObject, eventdata, handles, 1, 'sensor1Check')
% --

%% --- Executes on button press in sensor2Check.
function sensor2Check_Callback(hObject, eventdata, handles)
setVisSensors(hObject, eventdata, handles, 2, 'sensor2Check')
% --

%% --- Executes on button press in sensor3Check.
function sensor3Check_Callback(hObject, eventdata, handles)
setVisSensors(hObject, eventdata, handles, 3, 'sensor3Check')
% --

%% --- Executes on button press in sensor4Check.
function sensor4Check_Callback(hObject, eventdata, handles)
setVisSensors(hObject, eventdata, handles, 4, 'sensor4Check')
% --

%% Executes on checkbox select show load of Sensor#. Plots sensor force data in addition to the force data in Ax2
function setVisSensors(hObject, eventdata, handles, SensorNum, SensorName)
% -- retrieve data
h1 = handles.currentData.plotHandles.hLsensor(SensorNum,1);
h2 = handles.currentData.plotHandles.hLsensor(SensorNum,2);
h3 = handles.currentData.plotHandles.hLsensor(SensorNum,3);
% --
vis = 'off'; 
if get(handles.(SensorName),'Value'), vis = 'on'; end;
set(h1,'Visible',vis);
set(h2,'Visible',vis);
set(h3,'Visible',vis);
% --
vis = 'off';
checkSum = [ get(handles.sensor1Check,'Value');
             get(handles.sensor2Check,'Value');
             get(handles.sensor3Check,'Value');
             get(handles.sensor4Check,'Value'); ];
if sum( checkSum )>0
    vis = 'on';
end;
set(handles.currentData.plotHandles.panelPosHdl,'Visible',vis);
% -- set back data
handles.currentData.plotHandles.hLsensor(SensorNum,:) = [h1 h2 h3];
guidata(hObject,handles);
% -- End of setVisSensors

%% --- Executes on button press in forceOnlyCheck.
function forceOnlyCheck_Callback(hObject, eventdata, handles)
child3 = get(handles.ax3,'Children');
ax3R    = findobj('Tag','ax3R');
child3R = get(ax3R,'Children');
child4 = get(handles.ax4,'Children');
vis = 'on';
if get(handles.forceOnlyCheck,'Value'), vis = 'off'; end;
try set(handles.ax3,'Visible',vis); end; %#ok<*TRYNC>
for i=1:numel(child3)
    try set(child3(i),'Visible',vis); end;
end;
try set(ax3R,'Visible',vis); end; %#ok<*TRYNC>
for i=1:numel(child3R)
    try set(child3R(i),'Visible',vis); end;
end;
try  set(handles.ax4,'Visible',vis); end;
for i=1:numel(child4)
    try set(child4(i),'Visible',vis); end;
end;
if ~get(handles.forceOnlyCheck,'Value')
    h = findobj('Tag','joiaTimeTraceCheck');
    callback = get(h,'Callback');
        feval(callback,hObject,eventdata);
    h = findobj('Tag','sensor1Check');
    callback = get(h,'Callback');
        feval(callback,hObject,eventdata);
    h = findobj('Tag','sensor2Check');
    callback = get(h,'Callback'); 
        feval(callback,hObject,eventdata);
    h = findobj('Tag','sensor3Check');
    callback = get(h,'Callback');
        feval(callback,hObject,eventdata);
    h = findobj('Tag','sensor4Check');
    callback = get(h,'Callback');
        feval(callback,hObject,eventdata);   
    set(handles.ax2,'Position', [handles.currentData.plotHandles.origPosAx2],'XTickLabel','');
else
    ax2Top      = handles.currentData.plotHandles.origPosAx2(2) + handles.currentData.plotHandles.origPosAx2(4);
    ax4Bot      = handles.currentData.plotHandles.origPosAx4(2);
    newHeight   = ax2Top - ax4Bot;
    sameWidth   = handles.currentData.plotHandles.origPosAx2(3);
    sameLeft    = handles.currentData.plotHandles.origPosAx2(1);
    set(handles.ax2,'Position', [sameLeft  ax4Bot  sameWidth  newHeight],'XTickLabel',get(handles.ax4,'XTickLabel'));
end;
% -- End of forceOnlyCheck_Callback

%% --- Executes on button press in useSmoothData.
function useSmoothData_Callback(hObject, eventdata, handles)

%%
function avgWindowSize_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function avgWindowSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
function EdStart_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function EdStart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
function EdEnd_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function EdEnd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
function saveHPZstats_Callback(hObject, eventdata, handles)

Data            = handles.currentData;
DeltaTime       = Data.TactileSensor.ASFfile.SECONDS_PER_FRAME;
X               = Data.calculations.X;
Y               = Data.calculations.Y;
Pressures       = Data.TactileSensor.ASFfile.Pressures;
timeVec         = Data.calculations.timeVec;
currentFolder   = handles.options.currentFolder;
% --
selectedObj = get(handles.saveHPZopts,'SelectedObject');
if selectedObj==handles.saveAllOpt
    iLim = [1 numel(timeVec)];
elseif selectedObj==handles.saveRangeOpt
    iLim = [ get(handles.EdStart,'UserData')  get(handles.EdEnd,'UserData') ];
end
% --
numFramesTotal  = size(Pressures,3);
numFrames2do    = (iLim(2)-iLim(1)+1);
blankCell       = cell(numFramesTotal,1);
% --
hWait = waitbar(0,'Calculations on each frame in progress ...');
for i = iLim(1):iLim(2)
    % --
    waitbar( (i-iLim(1))/numFrames2do )
    % -- data from each frame
    foundHPZs_data = findHPZs( X,Y,Pressures(:,:,i),handles );
    if i==iLim(1)
        % initialize the output structure dynamically with all the content of foundHPZs_data
        FieldNames = fieldnames(foundHPZs_data);
        for iStruct = 1:numel(FieldNames)
            if strcmpi(FieldNames{iStruct},'correlationData')
                % --
                blank_thetaCell   	= cell( size(foundHPZs_data.bins_theta,1),1 );
                correlationData 	= blank_thetaCell;
                blank_distCell      = cell( size(foundHPZs_data.bins_D,1),1 );
                null_distMatrix     = zeros( size(foundHPZs_data.bins_D,1),1 );
                for iCorr = 1:numel(correlationData)
                    correlationData{iCorr} = blank_distCell;
                    for jCorr = 1:size(foundHPZs_data.bins_D,1)
                        correlationData{iCorr}{jCorr,1} = blank_thetaCell;
                        for kCorr = 1:size(blank_thetaCell,1)
                            correlationData{iCorr}{jCorr,1}{kCorr,1} = null_distMatrix;
                        end
                    end
                end
                % --
                calculatedData.(FieldNames{iStruct}) = correlationData;
            else
                calculatedData.(FieldNames{iStruct}) = blankCell;
            end
        end
    end
    % -- fill in new values in output structure
    for iStruct = 1:numel(FieldNames)
        if strcmpi(FieldNames{iStruct},'correlationData')
            % --
            existingCorrelationData = calculatedData.(FieldNames{iStruct});
            correlationData = foundHPZs_data.(FieldNames{iStruct});
            for iCorr = 1:numel(correlationData)
                for jCorr = 1:numel(correlationData{1})
                    for kCorr = 1:numel(correlationData{2})
                     	correlationData{iCorr}{jCorr,1}{kCorr,1} = existingCorrelationData{iCorr}{jCorr,1}{kCorr,1} + correlationData{iCorr}{jCorr,1}{kCorr,1};
                    end
                end
            end
            % --
            calculatedData.(FieldNames{iStruct}) = correlationData;
        elseif ismember(FieldNames{iStruct},{'bins_theta','bins_D'})
            if i==iLim(1)
                calculatedData.(FieldNames{iStruct}) = foundHPZs_data.(FieldNames{iStruct});
            end
        else
            calculatedData.(FieldNames{iStruct}){i} = foundHPZs_data.(FieldNames{iStruct});
        end
    end
    % --
end
close(hWait); pause(0.1);
% -- get time traces for individual HPZ
individualHPZs_data = getTimeTracesHPZs( calculatedData,X,Y,DeltaTime,true ); %#ok<NASGU>
% -- 
% save but first, look if file already exists
listExistingTemp = dir('*.mat'); listExisting = cell(numel(listExistingTemp),1);
for i=1:numel(listExistingTemp), listExisting{i} = listExistingTemp(i).name; end
currentName = ['calculatedData_',currentFolder,'.mat'];
ok = 0; count = 1;
while ~ok
    if ismember(currentName,listExisting)
        if count==1, currentName = [currentName(1:end-4),' (',num2str(count),')',currentName(end-3:end)];
        else
            ind(1) = strfind(currentName,'(');   ind(2) = strfind(currentName,')');
            currentName = [currentName(1:ind(1)),num2str(count),currentName(ind(2):end)];
        end
    else ok = 1;
    end
    count = count+1;
end
save( currentName,'calculatedData','individualHPZs_data' );
% -- turn objects enabling on
SetObjOnOff(handles,'on')
SetObjOnOff(handles,'on',{'stopButton'})
% -- End of SaveHPZ


%% --- Executes when selected object is changed in saveHPZopts.
function saveHPZopts_SelectionChangeFcn(hObject, eventdata, handles)
% --


%% Defaults Callback functions for MATLAB
function mult_std_peaks_Callback(hObject, eventdata, handles)
EnableAx10_Callback(hObject, eventdata, handles)
% --

%%
function mult_std_limits_Callback(hObject, eventdata, handles)
EnableAx10_Callback(hObject, eventdata, handles)
% --

%%
function minValue_peaks_Callback(hObject, eventdata, handles)
EnableAx10_Callback(hObject, eventdata, handles)
% --

function edit19_Callback(hObject, eventdata, handles)
EnableAx10_Callback(hObject, eventdata, handles)
% -- End of defalut callback functions


