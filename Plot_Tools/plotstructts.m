function varargout = plotstructts(varargin)
% PLOTSTRUCTTS MATLAB code for plotstructts.fig
%      PLOTSTRUCTTS, by itself, creates a new PLOTSTRUCTTS or raises the existing
%      singleton*.
%
%      H = PLOTSTRUCTTS returns the handle to a new PLOTSTRUCTTS or the handle to
%      the existing singleton*.
%
%      PLOTSTRUCTTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTSTRUCTTS.M with the given input arguments.
%
%      PLOTSTRUCTTS('Property','Value',...) creates a new PLOTSTRUCTTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotstructts_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotstructts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotstructts

% Last Modified by GUIDE v2.5 01-Dec-2015 16:43:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotstructts_OpeningFcn, ...
                   'gui_OutputFcn',  @plotstructts_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before plotstructts is made visible.
function plotstructts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotstructts (see VARARGIN)

% Choose default command line output for plotstructts
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using plotstructts.

if strcmp(get(hObject,'Visible'),'off')
    axis off
end
set(handles.axes1, 'Visible', 'off')
for i = 1:length(varargin)
    % 1. Safe data in handles
    handles.mydata{i} = varargin{i};     
end

handles.nr_data = i;
Variables = fieldnames(varargin{1}.Variables);

for i = 1:length(Variables)
    isfix(i) = isfixedvar(Variables{i});
end
Variables(isfix == 1) = [];

set(handles.Variable_selection, 'String', Variables);

if isfield(varargin{1}.Data, 'region_names')
    set(handles.Region_selection, 'String', varargin{1}.Data.region_names);
elseif isfield(varargin{1}.Data, 'station_name')
    set(handles.Region_selection, 'String', varargin{1}.Data.station_name);
elseif isfield(varargin{1}.Data, 'regions')
    set(handles.Region_selection, 'String', num2str(varargin{1}.Data.regions));
else
    set(handles.Region_selection, 'String', 'Region_1')
end

if isfield(varargin{1}.Data, 'regions')
    set(handles.Region_ID, 'String', num2str(varargin{1}.Data.regions));
    
   % tmp = mat2cell(1:length(varargin{1}.Data.regions), length(varargin{1}.Data.regions), 1);
    
elseif isfield(varargin{1}.Data, 'stations')
    set(handles.Region_ID, 'String', num2str(varargin{1}.Data.stations));
end




load coast

handles.coast.lat  = lat;
handles.coast.long = long;

axes(handles.Station_Map);

plot(long, lat, 'k', 'linewidth', 1.5)
set(handles.Station_Map, 'Ytick', []);
set(handles.Station_Map, 'Xtick', []);



 guidata(hObject,handles);                  
                    
% UIWAIT makes plotstructts wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plotstructts_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


axes(handles.axes1);

set(handles.axes1, 'Visible', 'off')

cla;


% Get user parameter

% List of variables in the input data
variables = get(handles.Variable_selection, 'String');  
% ID of the selected variable
var_id    = get(handles.Variable_selection, 'Value');
% Name of the selected variable
plotvar   = variables(var_id);

% List of regions in the input data
regions   = get(handles.Region_selection, 'String');  
% ID of the selected region
reg_id    = get(handles.Region_selection, 'Value');  
% Name of the selected region
plotreg   = regions(reg_id);

% Y-limits
y_min = str2num(get(handles.Y_min, 'String'));
y_max = str2num(get(handles.Y_max, 'String'));

% 
stats = get(handles.Stat_box, 'Value');

% Time limits
year_min  = str2num(get(handles.Year_min, 'String'));
month_min = str2num(get(handles.Month_min, 'String'));
day_min   = str2num(get(handles.Day_min, 'String'));

year_max  = str2num(get(handles.Year_max, 'String'));
month_max = str2num(get(handles.Month_max, 'String'));
day_max   = str2num(get(handles.Day_max, 'String'));

if isempty(y_min)
    if find(ismember(handles.mydata{1}.Variables.(char(plotvar)).dimensions, ...
            'time')) == 1
        y_min = min(handles.mydata{1}.Data.(char(plotvar))(:, reg_id));
    elseif find(ismember(handles.mydata{1}.Variables.(char(plotvar)).dimensions, ...
            'time')) == 2
        y_min = min(handles.mydata{1}.Data.(char(plotvar))(reg_id, :));
    end
end

if isempty(y_max) 
     if find(ismember(handles.mydata{1}.Variables.(char(plotvar)).dimensions, ...
            'time')) == 1
        y_max = max(handles.mydata{1}.Data.(char(plotvar))(:, reg_id));
    elseif find(ismember(handles.mydata{1}.Variables.(char(plotvar)).dimensions, ...
            'time')) == 2
        y_max = max(handles.mydata{1}.Data.(char(plotvar))(reg_id, :));
    end
end

if isempty(year_min)
    tme_start = datevec(handles.mydata{1}.TimeStamp(1, :));
    year_min  = tme_start(1);
end

if isempty(month_min), month_min = 1; end;
if isempty(day_min), day_min = 1; end

if isempty(year_max)
    tme_end = datevec(handles.mydata{1}.TimeStamp(end, :));
    year_max = tme_end(1);
end

if isempty(month_max), month_max = 12; end;
if isempty(day_max), day_max = eomday(year_max, month_max); end

lnewdth = str2num(get(handles.Linewidth, 'String'));
if isempty(lnewdth)
    lnewdth = 1.5;
end

% First, extract the map at the specified date from the data
for i = 1:handles.nr_data

%     if stats == 1  
%         if i > 1
%             handles.mydata{i} = trunc_TS(handles.mydata{i}, ...
%                                          handles.mydata{1}.Data.time(1, :), ...
%                                          handles.mydata{1}.Data.time(end, :));
%         end
%     end           
    if isfield(handles.mydata{i}.DataInfo, 'title')
        ttle{i} = handles.mydata{i}.DataInfo.title;
    else
        ttle{i} = ['Dataset ', num2str(i)];
    end
    
    tme_indx = ...
        find(ismember(handles.mydata{i}.Variables.(char(plotvar)).dimensions, 'time') ...
                                                                     == 1);
                                                                 
    if tme_indx == 1                                                             
        plotdata{i} = handles.mydata{i}.Data.(char(plotvar))(:, reg_id);
    elseif tme_indx == 2
        plotdata{i} = handles.mydata{i}.Data.(char(plotvar))(reg_id, :);
    end
    
    hc = plot(datetime(handles.mydata{i}.Data.time), plotdata{i}, 'linewidth', lnewdth);
        
    hold on
end

if ~isempty(handles.mydata{1}.Variables.(char(plotvar)).units)
    set(get(handles.axes1, 'Ylabel'), 'String', handles.mydata{1}.Variables.(char(plotvar)).units, 'fontsize', 16);
end

ylim([y_min y_max]);
T_min = datenum(year_min, month_min, day_min);
T_max = datenum(year_max, month_max, day_max);
xlim([T_min T_max])
%set(get(handles.axes1, 'YLim'), [y_min y_max])


% % Truncate the spatial extent map according to the data
% axis([min(lons) max(lons) min(lats) max(lats)]) 
lgnd = legend(ttle);
set(lgnd, 'Fontsize', 14);
% Hold off for new maps!
hold off







% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)








% --- Executes on selection change in Variable_selection.
function Variable_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Variable_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Variable_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Variable_selection


% --- Executes during object creation, after setting all properties.
function Variable_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Variable_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Region_selection.
function Region_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Region_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Region_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Region_selection

reg_id    = get(handles.Region_selection, 'Value');  
set(handles.Region_ID, 'Value', reg_id);

axes(handles.Station_Map);
plot(handles.coast.long, handles.coast.lat, 'k', 'linewidth', 1.5)
set(handles.Station_Map, 'Ytick', []);
set(handles.Station_Map, 'Xtick', []);
hold on
if isfield(handles.mydata{1}.Data, 'lon') & ...
                                     isfield(handles.mydata{1}.Data, 'lat') 
    plot(handles.mydata{1}.Data.lon(reg_id), ...
                handles.mydata{1}.Data.lat(reg_id), 'r.', 'markersize', 20)
end
hold off
guidata(hObject,handles);    


% --- Executes during object creation, after setting all properties.
function Region_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Region_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_min_Callback(hObject, eventdata, handles)
% hObject    handle to Y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y_min as text
%        str2double(get(hObject,'String')) returns contents of Y_min as a double


% --- Executes during object creation, after setting all properties.
function Y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_max_Callback(hObject, eventdata, handles)
% hObject    handle to Y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y_max as text
%        str2double(get(hObject,'String')) returns contents of Y_max as a double


% --- Executes during object creation, after setting all properties.
function Y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Month_min_Callback(hObject, eventdata, handles)
% hObject    handle to Month_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Month_min as text
%        str2double(get(hObject,'String')) returns contents of Month_min as a double


% --- Executes during object creation, after setting all properties.
function Month_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Month_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Day_min_Callback(hObject, eventdata, handles)
% hObject    handle to Day_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Day_min as text
%        str2double(get(hObject,'String')) returns contents of Day_min as a double


% --- Executes during object creation, after setting all properties.
function Day_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Day_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Year_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function Year_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Year_max_Callback(hObject, eventdata, handles)
% hObject    handle to Year_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Year_max as text
%        str2double(get(hObject,'String')) returns contents of Year_max as a double


% --- Executes during object creation, after setting all properties.
function Year_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Year_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Month_max_Callback(hObject, eventdata, handles)
% hObject    handle to Month_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Month_max as text
%        str2double(get(hObject,'String')) returns contents of Month_max as a double


% --- Executes during object creation, after setting all properties.
function Month_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Month_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Day_max_Callback(hObject, eventdata, handles)
% hObject    handle to Day_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Day_max as text
%        str2double(get(hObject,'String')) returns contents of Day_max as a double


% --- Executes during object creation, after setting all properties.
function Day_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Day_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Linewidth_Callback(hObject, eventdata, handles)
% hObject    handle to Linewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Linewidth as text
%        str2double(get(hObject,'String')) returns contents of Linewidth as a double


% --- Executes during object creation, after setting all properties.
function Linewidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Linewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Region_ID.
function Region_ID_Callback(hObject, eventdata, handles)
% hObject    handle to Region_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Region_ID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Region_ID

reg_id    = get(handles.Region_ID, 'Value')  

set(handles.Region_selection, 'Value', reg_id);

axes(handles.Station_Map);
plot(handles.coast.long, handles.coast.lat, 'k', 'linewidth', 1.5)
set(handles.Station_Map, 'Ytick', []);
set(handles.Station_Map, 'Xtick', []);
hold on
if isfield(handles.mydata{1}.Data, 'lon') & ...
                                     isfield(handles.mydata{1}.Data, 'lat') 
    plot(handles.mydata{1}.Data.lon(reg_id), ...
                handles.mydata{1}.Data.lat(reg_id), 'r.', 'markersize', 20)
end
hold off

guidata(hObject,handles);    


% --- Executes during object creation, after setting all properties.
function Region_ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Region_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in Stat_box.
function Stat_box_Callback(hObject, eventdata, handles)
% hObject    handle to Stat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Stat_box
