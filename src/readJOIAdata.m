function Data = read_JOIAdata(Folder)
%READJOIADATA Reads JOIA data.
%   DATA = READJOIADATA('C:\theFolder') reads JOIA data available in 
%   'C:\theFolder' and returns it in the structure named DATA
%
%   JOIA data in folders normally includes following files:
%   -       *.asf: Tactile sensors data
%   -  *_area.asg: Areas calculated from tactile sensors data
%   - *_force.asg: Forces calculated from tactile sensors data
%   -       *.cal: Calibration
%   -  *_(kN).txt: Forces from load cells
%   - *_(MPa).txt: Pressures from load cells
%   -       *.fbx: Proprietary format file from TEKSCAN (Not read)
%   -       *.FSX: Proprietary format file from TEKSCAN(Not read)
% 
%   If one or more file (as listed above) are missing, the output for that
%   file will be an empty value.
% 
%   the output variable DATA is a structure containg following main fields:
%   - DATA.TACTILESENSOR.ASFFILE: Contains data from *.asf file.
%   - DATA.TACTILESENSOR.ASGFILE: Contains data from both *_area.asg and *_force.asg files.
%   - DATA.TACTILESENSOR.CALFILE: Contains data from *.txt file.
%   - DATA.LOADCELLS.TXTFILE: Contains data from both *_(kN).txt and *_(MPa).txt files.
% 
%   EXAMPLE:
%   Data = readJOIAdata( 'C:\Users\mrichard\Documents\270851-RDC IGNITE High Pressure Zones\Local\JOIA\Data\Data_1998\980126-1' );
% 
% 
% 
%   Martin Richard
%   C-CORE, St. John's, NL CANADA


%% Allowed file types and init.
AllowedTypes = {'asf','_area.asg','_force.asg','cal','_(kN).txt','_(MPa).txt'};
Data = struct('TactileSensor',[], 'LoadCells',[]);
    Data.TactileSensor.ASFfile = [];
    Data.TactileSensor.ASGfile = [];
    Data.TactileSensor.CALfile = [];
    Data.LoadCells.TXTfile     = [];
if strcmp(Folder(end),filesep), Folder = Folder(1:end-1); end


%% Read sequentially
for i=1:numel(AllowedTypes)
    gotErrors = false;
    CurrentFile = dir( [Folder,filesep,'*',AllowedTypes{i}] );
    if ~isempty(CurrentFile) && numel(CurrentFile)==1
        switch AllowedTypes{i}
            case 'asf'          , [Data gotErrors] = readASF( Data,[Folder,filesep,CurrentFile.name] );
            case '_area.asg'    , [Data gotErrors] = readASG( Data,[Folder,filesep,CurrentFile.name],'area' );
            case '_force.asg'   , [Data gotErrors] = readASG( Data,[Folder,filesep,CurrentFile.name],'force' );
            case 'cal'          , [Data gotErrors] = readCAL( Data,[Folder,filesep,CurrentFile.name] );
            case '_(kN).txt'    , [Data gotErrors] = readTXT( Data,[Folder,filesep,CurrentFile.name],'force' );
            case '_(MPa).txt'   , [Data gotErrors] = readTXT( Data,[Folder,filesep,CurrentFile.name],'press' );
        end
    else
        gotErrors = true;
    end
    if gotErrors
        mssg = 'File is either missing or there was more than one file of this type';
        switch AllowedTypes{i}
            case 'asf',         Data.TactileSensor.ASFfile = mssg;
            case '_area.asg',   Data.TactileSensor.ASGfile.area = mssg;
            case '_force.asg', 	Data.TactileSensor.ASGfile.force = mssg;
            case 'cal',         Data.TactileSensor.CALfile = mssg;
            case '_(kN).txt',   Data.LoadCells.TXTfile.force = mssg;
            case '_(MPa).txt',	Data.LoadCells.pressure = mssg;
        end
    end
end


%% sub-function readASF
function [Data gotErrors] = readASF( Data,FileName )
% --
[HeaderData gotErrors] = readGenHeader( FileName );
if gotErrors, return; end;
Data.TactileSensor.ASFfile = HeaderData;
% --
fid = fopen(FileName);
if fid<0, gotErrors = true; return; end;
% --
isHeader = 1;
tline = fgetl(fid);
while ischar(tline)
    if isempty(tline), isHeader = 0; end
    % --
    tline = fgetl(fid);
    % --
    if ~isHeader
        % --
        str4reading = '';
        for iColumns = 1:Data.TactileSensor.ASFfile.COLS
            str4reading = [str4reading, ' %f'];
        end
        % --
        TactileData = nan( Data.TactileSensor.ASFfile.ROWS , Data.TactileSensor.ASFfile.COLS , Data.TactileSensor.ASFfile.END_FRAME-Data.TactileSensor.ASFfile.START_FRAME+1 );
        % --
        hWait = waitbar(0,'Reading ASF file, please wait ...');
        for iFrame = Data.TactileSensor.ASFfile.START_FRAME:Data.TactileSensor.ASFfile.END_FRAME
            waitbar(iFrame/(Data.TactileSensor.ASFfile.END_FRAME-Data.TactileSensor.ASFfile.START_FRAME+1))
            CurrentFrame = textscan(fid, str4reading , Data.TactileSensor.ASFfile.ROWS , 'delimiter',',' , 'CollectOutput',1 );
            TactileData(:,:,iFrame-Data.TactileSensor.ASFfile.START_FRAME+1) = flipud(CurrentFrame{1});
            tline = fgetl(fid);
            while isempty(strfind(tline,'Frame'))
                tline = fgetl(fid);
                if ~isempty(strfind(tline,Data.TactileSensor.ASFfile.ASCII_DATA))
                    tline = []; %  to break from main loop.
                    break
                end
            end
        end
        % --
        close(hWait)
        Data.TactileSensor.ASFfile.Pressures = TactileData;
    end
end
fclose(fid);


%% sub-function readASG
function [Data gotErrors] = readASG( Data,FileName,FileType )
% --
[HeaderData gotErrors] = readGenHeader( FileName );
if gotErrors, return; end;
Data.TactileSensor.ASGfile.(FileType) = HeaderData;
% --
fid = fopen(FileName);
if fid<0, gotErrors = true; return; end;
% --
isHeader = 1;
tline = fgetl(fid);
while ischar(tline)
    if strcmpi(tline,'ASCII_DATA @@'), isHeader = 0; end
    % --
    tline = fgetl(fid);
    % --
    if ~isHeader
        % --
        ind = strfind(tline,', ');
        Fields{1,1} = tline(1:ind(1)-1); %#ok<*AGROW>
        for i=2:numel(ind)
            Fields{i,1} = tline(ind(i-1)+2:ind(i)-1);
        end
        Fields{end+1,1} = tline(ind(end)+2:end);
        % --
        str4reading = '';
        for iColumns = 1:numel(Fields)
            str4reading = [str4reading, ' %f']; %#ok<AGROW>
        end
        % --
        CurrentFrame = textscan(fid, str4reading , Data.TactileSensor.ASFfile.END_FRAME-Data.TactileSensor.ASFfile.START_FRAME+1 , 'delimiter',',' , 'CollectOutput',1 );
        for i=1:numel(Fields)
            ind = strfind(Fields{i},'(');
            Data.TactileSensor.ASGfile.(FileType).(Fields{i}(1:ind-1)) = CurrentFrame{1}(:,i);
        end
        tline = []; % to break the while loop
    end
end
fclose(fid);


%% sub-function readCAL
function [Data gotErrors] = readCAL( Data,FileName )
% --
[HeaderData gotErrors] = readGenHeader( FileName );
if gotErrors, return; end;
Data.TactileSensor.CALfile = HeaderData;


%% sub-function readTXT
function [Data gotErrors] = readTXT( Data,FileName,FileType )
% --
gotErrors = false;
% --
switch FileType
    case 'force'
        Headers = {'PlottingFreq','ElapsedTime_sec','PenetrationQuantity_cm','Fz1_kN','Fz2_kN','Fz3_kN','Fz4_kN','Fz5_kN','Fz6_kN','Fz7_kN',...
            'Fz8_kN','Fz9_kN','Fz10_kN','Fz11_kN','Fz12_kN','Fz13_kN','Fz14_kN','Fz15_kN','Anglemeter_deg','Accelorometer_gal',...
            'Displacement_mm','EntireLoadA_kN','Lz4All_kN','ActionPoint_cm','HydraulicMeasFixedCommand_kgfPERcmPERcm','EntireLoadOilPressGauge_kN',...
            'L4Fx4_kN','L4Fx5_kN','L4Fx6_kN','L4Fx7_kN','L4Fx8_kN','L4Fx9_kN','L4Fx10_kN','L4Fx11_kN','L4Fx12_kN','L4Fx13_kN',...
            'L4Fy4_kN','L4Fy5kN','L4Fy6_kN','L4Fy7_kN','L4Fy8_kN','L4Fy9_kN','L4Fy10_kN','L4Fy11_kN','L4Fy12_kN','L4Fy13_kN',...
            'L4My4_kNm','L4My5_kNm','L4My6_kNm','L4My7_kNm','L4My8_kNm','L4My9_kNm','L4My10_kNm','L4My11_kNm','L4My12_kNm','L4My13_kNm','BE'};
    case 'press'
        Headers = {'PlottingFreq','ElapsedTime_sec','PenetrationQuantity_cm','Fz1_MPa','Fz2_MPa','Fz3_MPa','Fz4_MPa','Fz5_MPa','Fz6_MPa','Fz7_MPa',...
            'Fz8_MPa','Fz9_MPa','Fz10_MPa','Fz11_MPa','Fz12_MPa','Fz13_MPa','Fz14_MPa','Fz15_MPa','Anglemeter_deg','Accelorometer_gal',...
            'Displacement_mm','EntireLoadA_MPa','Lz4All_MPa','ActionPoint_cm','HydraulicMeasFixedCommand_kgfPERcmPERcm','EntireLoadOilPressGauge_MPa',...
            'L4Fx4_MPa','L4Fx5_MPa','L4Fx6_MPa','L4Fx7_MPa','L4Fx8_MPa','L4Fx9_MPa','L4Fx10_MPa','L4Fx11_MPa','L4Fx12_MPa','L4Fx13_MPa',...
            'L4Fy4_MPa','L4Fy5MPa','L4Fy6_MPa','L4Fy7_MPa','L4Fy8_MPa','L4Fy9_MPa','L4Fy10_MPa','L4Fy11_MPa','L4Fy12_MPa','L4Fy13_MPa',...
            'L4My4_kNm','L4My5_kNm','L4My6_kNm','L4My7_kNm','L4My8_kNm','L4My9_kNm','L4My10_kNm','L4My11_kNm','L4My12_kNm','L4My13_kNm','BE'};
end
% --
Matrix = dlmread(FileName, '\t', 1, 0);
if size(Matrix,2)>numel(Headers), tooBig = size(Matrix,2)-numel(Headers) - 1; Matrix(:,end-tooBig:end) = []; end
% --
for i=1:size(Matrix,2)
    Data.LoadCells.TXTfile.(FileType).(Headers{i}) = Matrix(:,i);
end


%% sub-function readGenHeader
function [HeaderData gotErrors] = readGenHeader( FileName )
gotErrors = false;
% --
fid = fopen(FileName);
if fid<0, gotErrors = true; return; end;
% --
isHeader = 1;
while isHeader
    tline = fgetl(fid);
    % --
    if strcmpi(tline,'ASCII_DATA @@') || isempty(tline) || ~ischar(tline), isHeader = 0; end
    % --
    if ischar(tline)
        if ~isempty(strfind(tline,'First name:'))
            ind = strfind(tline,':'); ind = ind(1);
            HeaderData.First_name = tline(ind+2:end);
        elseif ~isempty(strfind(tline,'Last name:'))
            ind = strfind(tline,':'); ind = ind(1);
            HeaderData.Last_name = tline(ind+2:end);
        elseif ~isempty(strfind(tline,'Date:'))
            ind = strfind(tline,':'); ind = ind(1);
            HeaderData.Date = tline(ind+2:end);
            ind = [strfind(HeaderData.Date,'/') strfind(HeaderData.Date,' ') strfind(HeaderData.Date,':')];
            mm = str2num(HeaderData.Date(1:ind(1)-1)); %#ok<*ST2NM>
            dd = str2num(HeaderData.Date(ind(1)+1:ind(2)-1));
            yy = str2num(['19',HeaderData.Date(ind(2)+1:ind(3)-1)]);
            hh = str2num(HeaderData.Date(ind(3)+1:ind(4)-1));
            min= str2num(HeaderData.Date(ind(4)+1:end));
            HeaderData.Date_MATLAB = datenum(yy,mm,dd,hh,min,0);
        else
            ind = strfind(tline,' ');
            if ismember(tline(1:ind(1)-1),{'ROWS','COLS','NOISE_THRESHOLD','SECONDS_PER_FRAME','MICRO_SECOND','START_FRAME','END_FRAME'})
                HeaderData.(tline(1:ind(1)-1)) = str2num( tline(ind(1)+1:end) );
            elseif ismember(tline(1:ind(1)-1),{'ROW_SPACING','COL_SPACING','SENSEL_AREA','SATURATION_PRESSURE'}) && ~strcmpi(FileName(end-2:end),'cal')
                HeaderData.(tline(1:ind(1)-1)) = str2num( tline(ind(1)+1:ind(2)) );
                HeaderData.([tline(1:ind(1)-1),'_units']) = tline(ind(2)+1:end);
            elseif ismember(tline(1:ind(1)-1),{'CALIBRATION_POINT_1'})
                HeaderData.(tline(1:ind(1)-1)).Force              = str2num( tline(ind(1)+1:ind(2)) );
                HeaderData.(tline(1:ind(1)-1)).Force_units        = tline(ind(2)+2:ind(3)-2);
                HeaderData.(tline(1:ind(1)-1)).RawSum             = str2num( tline(ind(3)+1:ind(4)) );
                HeaderData.(tline(1:ind(1)-1)).NumberLoadedCells  = str2num( tline(ind(6)+1:ind(7)) );
            else
                HeaderData.(tline(1:ind(1)-1)) = tline(ind(1)+1:end);
            end
        end
    end
    % --
end
fclose(fid);
