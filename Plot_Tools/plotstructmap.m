function varargout = plotstructmap(varargin)
% PLOTSTRUCTMAP MATLAB code for plotstructmap.fig
%      PLOTSTRUCTMAP, by itself, creates a new PLOTSTRUCTMAP or raises the existing
%      singleton*.
%
%      H = PLOTSTRUCTMAP returns the handle to a new PLOTSTRUCTMAP or the handle to
%      the existing singleton*.
%
%      PLOTSTRUCTMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTSTRUCTMAP.M with the given input arguments.
%
%      PLOTSTRUCTMAP('Property','Value',...) creates a new PLOTSTRUCTMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotstructmap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotstructmap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotstructmap

% Last Modified by GUIDE v2.5 30-Oct-2015 18:00:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotstructmap_OpeningFcn, ...
                   'gui_OutputFcn',  @plotstructmap_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
schnupfen
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before plotstructmap is made visible.
function plotstructmap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotstructmap (see VARARGIN)

% Choose default command line output for plotstructmap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using plotstructmap.

if strcmp(get(hObject,'Visible'),'off')
    axis off
end

handles.mydata = varargin{1};

if isfield(varargin{1}.Variables, 'time')
    datestrings = datestr(varargin{1}.Data.time);
else
    datestrings = {'None'};
end

Variables      = fieldnames(varargin{1}.Variables);

for i = 1:length(Variables)
    isfix(i) = isfixedvar(Variables{i});
    
    if strcmp(Variables{i}, 'region_map')
        isfix(i) = 0;
    end
end
Variables(isfix == 1) = [];



set(handles.popupmenu1, 'String', datestrings);
set(handles.popupmenu2, 'String', Variables);

set(handles.popupmenu3, 'String', ...
                        {'Imagesc', ...
                         'Pcolor', ...
                         'Contour', ...
                         'Contourf'});
                     
set(handles.popupmenu4, 'String', {'parula', ...
                        'jet', ...
                        'hsv', ...
                        'hot', ...
                        'cool', ...
                        'spring', ...
                        'summer', ...
                        'autumn', ...
                        'winter', ...
                        'gray', ...
                        'bone', ...
                        'copper', ...
                        'pink', ...
                        'lines', ...
                        'colorcube', ...
                        'prism', ...
                        'flag', ...
                        'white', ...
                        'redbluecmap'});

 guidata(hObject,handles);                  
                    
% UIWAIT makes plotstructmap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plotstructmap_OutputFcn(hObject, eventdata, handles)
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

cla;

% Get user parameter

% Index of the date which should be plotted


if strcmp(get(handles.popupmenu1, 'String'), 'none')
    date_index = 1;
else
   date_index = get(handles.popupmenu1, 'Value');   % Date
end

% List of variables in the input data
variables  = get(handles.popupmenu2, 'String');  
% ID of the selected variable
var_id     = get(handles.popupmenu2, 'Value');
% Name of the selected variable
plotvar    = variables(var_id);

% ID of the Plottype (imagesc, pcolor, contour, contourf)
plttype = get(handles.popupmenu3, 'Value');

% List of available colormaps
clrmps  = get(handles.popupmenu4, 'String');
% ID of the selected colormap
clr_id  = get(handles.popupmenu4, 'Value');
% Name of the selected colormap
clrmps  = clrmps{clr_id, :};


% First, extract the map at the specified date from the data
plotdata = handles.mydata.Data.(char(plotvar));
dims     = size(plotdata);

if length(dims) > 3
    % Data has different levels
    plotdata = squeeze(plotdata(date_index, 1, :, :));
    %levs = handles.mydata.(Data{2});
    
elseif length(dims) == 3
    % Data has 3 dimensions --> time, lat, lon
    plotdata = squeeze(plotdata(date_index, :, :));
end

% Get latitudes and longitudes of the data
lats = handles.mydata.Data.lat;
lons = handles.mydata.Data.lon;

% Should the function plot some boarder-lines as well?
cstlns   = get(handles.checkbox1, 'Value'); % Coastlines
brdrs    = get(handles.checkbox2, 'Value'); % National
ctchmnts = get(handles.checkbox3, 'Value'); % Catchments

 % Min/max value of the color axis
clr_min = str2num(get(handles.edit1, 'String'));
clr_max = str2num(get(handles.edit2, 'String'));

% If the boxes are empty, set the caxis to the max/min values in the data
    if isempty(clr_min) && isempty(clr_max) 
        cxis = [min(min(plotdata)) max(max(plotdata))];
    elseif clr_min == clr_max
        warning('off', 'backtrace')
        warning('Color data min and max must have different values!')
        cxis = [min(min(plotdata)) max(max(plotdata))];
    elseif clr_min > clr_max
        warning('off', 'backtrace')
        warning('Color data min must be smaller than color data max!')
        cxis = [min(min(plotdata)) max(max(plotdata))];
    else
        cxis = [clr_min, clr_max];
    end




% First, create a map with the specified method
switch plttype
	case 1
        hplot = imagesc(lons, lats, plotdata);
        set(hplot,'alphadata',~isnan(plotdata));
	case 2
        hplot = pcolor(lons, lats, plotdata);
        set(hplot,'alphadata',~isnan(plotdata));
        shading interp
	case 3 
        hplot = contour(lons, lats, plotdata);
    case 4
        hplot = contourf(lons, lats, plotdata);
end

% Hold plot for further options
hold on

% Add a colorbar to the map
hc = colorbar;

% Add some units to the colorbar
if isfield(handles.mydata.Variables.(char(plotvar)), 'units')
    set(get(hc, 'Ylabel'), 'String', handles.mydata.Variables.(char(plotvar)).units, 'fontsize', 16);
end


if cstlns == 1
    if ~isfield(handles, 'cstlns')
        load coast
        handles.cstlns = [long lat];
        % Update handles
        guidata(hObject, handles);
    end
	plot(handles.cstlns(:, 1), handles.cstlns(:, 2), 'k', 'linewidth', 2)
end
    
if ctchmnts == 1
    if ~isfield(handles, 'ctchlns')
        load ctchmnt_bnds_grdc.mat
        handles.ctchlns = ctchmnt_bnds;
        % Update handles
        guidata(hObject, handles);
    end
    plot(handles.ctchlns(:, 1), handles.ctchlns(:, 2), 'b', 'linewidth', 2)
end
    
if brdrs == 1
    if ~isfield(handles, 'brdrs')
        load natbrdrs_TMWB3.mat
        handles.brdrs = brdrs;
        % Update handles
        guidata(hObject, handles);
    end
    plot(handles.brdrs(:, 1), handles.brdrs(:, 2), 'r', 'linewidth', 2)
end

% Flip Y-Axis
axis xy


% Further options
% Set the selected colormap
colormap(clrmps) 
% Set the selected caxis
caxis(cxis)      
% Truncate the spatial extent map according to the data
axis([min(lons(:)) max(lons(:)) min(lats(:)) max(lats(:))]) 
  
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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);

% Check if we've reached the end of the timeseries
date_index = get(handles.popupmenu1, 'Value');

if date_index > 1
    cla;
    set(handles.popupmenu1, 'Value', get(handles.popupmenu1, 'Value') -1);

    % Update handles structure
    guidata(hObject, handles);

    % Get user parameter

    % Index of the date which should be plotted
    date_index = get(handles.popupmenu1, 'Value');   % Date

    % List of variables in the input data
    variables  = get(handles.popupmenu2, 'String');  
    % ID of the selected variable
    var_id     = get(handles.popupmenu2, 'Value');
    % Name of the selected variable
    plotvar    = variables(var_id);

    % ID of the Plottype (imagesc, pcolor, contour, contourf)
    plttype = get(handles.popupmenu3, 'Value');

    % List of available colormaps
    clrmps  = get(handles.popupmenu4, 'String');
    % ID of the selected colormap
    clr_id  = get(handles.popupmenu4, 'Value');
    % Name of the selected colormap
    clrmps  = clrmps{clr_id, :};


    % First, extract the map at the specified date from the data
    plotdata = handles.mydata.Data.(char(plotvar));
    dims     = size(plotdata);

    if length(dims) > 3
        % Data has different levels
        plotdata = squeeze(plotdata(date_index, 1, :, :));
        % Get latitudes and longitudes of the data
        lats = handles.mydata.Data.lat;
        lons = handles.mydata.Data.lon;
        levs = handles.mydata.(Data{2});
    
    elseif length(dims) == 3
        % Data has 3 dimensions --> time, lat, lon
        plotdata = squeeze(plotdata(date_index, :, :));
    
        % Get latitudes and longitudes of the data
        lats = handles.mydata.Data.lat;
        lons = handles.mydata.Data.lon;
    end


    % Should the function plot some boarder-lines as well?
    cstlns   = get(handles.checkbox1, 'Value'); % Coastlines
    brdrs    = get(handles.checkbox2, 'Value'); % National
    ctchmnts = get(handles.checkbox3, 'Value'); % Catchments

     % Min/max value of the color axis
    clr_min = str2num(get(handles.edit1, 'String'));
    clr_max = str2num(get(handles.edit2, 'String'));

    % If the boxes are empty, set the caxis to the max/min values in the data
    if isempty(clr_min) && isempty(clr_max) 
        cxis = [min(min(plotdata)) max(max(plotdata))];
    elseif clr_min == clr_max
        warning('off', 'backtrace')
        warning('Color data min and max must have different values!')
        cxis = [min(min(plotdata)) max(max(plotdata))];
    elseif clr_min > clr_max
        warning('off', 'backtrace')
        warning('Color data min must be smaller than color data max!')
        cxis = [min(min(plotdata)) max(max(plotdata))];
    else
        cxis = [clr_min, clr_max];
    end




    % First, create a map with the specified method
    switch plttype
        case 1
            hplot = imagesc(lons, lats, plotdata);
            set(hplot,'alphadata',~isnan(plotdata));
        case 2
            hplot = pcolor(lons, lats, plotdata);
            set(hplot,'alphadata',~isnan(plotdata));
            shading interp
        case 3 
            hplot = contour(lons, lats, plotdata);
        case 4
            hplot = contourf(lons, lats, plotdata);
    end

    % Hold plot for further options
    hold on

    % Add a colorbar to the map
    hc = colorbar;

    % Add some units to the colorbar
    if ~isempty(handles.mydata.Variables.(char(plotvar)).units)
        set(get(hc, 'Ylabel'), 'String', handles.mydata.Variables.(char(plotvar)).units, 'fontsize', 16);
    end


    if cstlns == 1
        if ~isfield(handles, 'cstlns')
            load coast
            handles.cstlns = [long lat];
            % Update handles
            guidata(hObject, handles);
        end
        plot(handles.cstlns(:, 1), handles.cstlns(:, 2), 'k', 'linewidth', 2)
    end
    
    if ctchmnts == 1
        if ~isfield(handles, 'ctchlns')
            load ctchmnt_bnds_grdc.mat
            handles.ctchlns = ctchmnt_bnds;
            % Update handles
            guidata(hObject, handles);
        end
        plot(handles.ctchlns(:, 1), handles.ctchlns(:, 2), 'b', 'linewidth', 2)
    end
    
    if brdrs == 1
        if ~isfield(handles, 'brdrs')
            load natbrdrs_TMWB3.mat
            handles.brdrs = brdrs;
            % Update handles
            guidata(hObject, handles);
        end
        plot(handles.brdrs(:, 1), handles.brdrs(:, 2), 'r', 'linewidth', 2)
    end
 
    % Flip Y-Axis
    axis xy

    % Further options
    % Set the selected colormap
    colormap(clrmps) 
    % Set the selected caxis
    caxis(cxis)      
    % Truncate the spatial extent map according to the data
    axis([min(lons(:)) max(lons(:)) min(lats(:)) max(lats(:))])
  
    % Hold off for new maps!
    hold off
    
else
    warning('off', 'backtrace')
    warning('plotstructmap.m: Reached first date of timeseries!')
end



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);

% Check if we've reached the end of the timeseries
date_index = get(handles.popupmenu1, 'Value');

if date_index < size(handles.mydata.Data.time, 1)
    cla;
    set(handles.popupmenu1, 'Value', get(handles.popupmenu1, 'Value') + 1);

    % Update handles structure
    guidata(hObject, handles);

    % Get user parameter

    % Index of the date which should be plotted
    date_index = get(handles.popupmenu1, 'Value');   % Date

    % List of variables in the input data
    variables  = get(handles.popupmenu2, 'String');  
    % ID of the selected variable
    var_id     = get(handles.popupmenu2, 'Value');
    % Name of the selected variable
    plotvar    = variables(var_id);

    % ID of the Plottype (imagesc, pcolor, contour, contourf)
    plttype = get(handles.popupmenu3, 'Value');

    % List of available colormaps
    clrmps  = get(handles.popupmenu4, 'String');
    % ID of the selected colormap
    clr_id  = get(handles.popupmenu4, 'Value');
    % Name of the selected colormap
    clrmps  = clrmps{clr_id, :};


    % First, extract the map at the specified date from the data
    plotdata = handles.mydata.Data.(char(plotvar));
    dims     = size(plotdata);

    if length(dims) > 3
        % Data has different levels
        plotdata = squeeze(plotdata(date_index, 1, :, :));
        % Get latitudes and longitudes of the data
        lats = handles.mydata.Data.lat;
        lons = handles.mydata.Data.lon;
        levs = handles.mydata.(Data{2});
    
    elseif length(dims) == 3
        % Data has 3 dimensions --> time, lat, lon
        plotdata = squeeze(plotdata(date_index, :, :));
    
        % Get latitudes and longitudes of the data
        lats = handles.mydata.Data.lat;
        lons = handles.mydata.Data.lon;
    end


    % Should the function plot some boarder-lines as well?
    cstlns   = get(handles.checkbox1, 'Value'); % Coastlines
    brdrs    = get(handles.checkbox2, 'Value'); % National
    ctchmnts = get(handles.checkbox3, 'Value'); % Catchments

     % Min/max value of the color axis
    clr_min = str2num(get(handles.edit1, 'String'));
    clr_max = str2num(get(handles.edit2, 'String'));

    % If the boxes are empty, set the caxis to the max/min values in the data
    if isempty(clr_min) && isempty(clr_max) 
        cxis = [min(min(plotdata)) max(max(plotdata))];
    elseif clr_min == clr_max
        warning('off', 'backtrace')
        warning('Color data min and max must have different values!')
        cxis = [min(min(plotdata)) max(max(plotdata))];
    elseif clr_min > clr_max
        warning('off', 'backtrace')
        warning('Color data min must be smaller than color data max!')
        cxis = [min(min(plotdata)) max(max(plotdata))];
    else
        cxis = [clr_min, clr_max];
    end




    % First, create a map with the specified method
    switch plttype
        case 1
            hplot = imagesc(lons, lats, plotdata);
            set(hplot,'alphadata',~isnan(plotdata));
        case 2
            hplot = pcolor(lons, lats, plotdata);
            set(hplot,'alphadata',~isnan(plotdata));
            shading interp
        case 3 
            hplot = contour(lons, lats, plotdata);
        case 4
            hplot = contourf(lons, lats, plotdata);
    end

    % Hold plot for further options
    hold on

    % Add a colorbar to the map
    hc = colorbar;

    % Add some units to the colorbar
    if ~isempty(handles.mydata.Variables.(char(plotvar)).units)
        set(get(hc, 'Ylabel'), 'String', handles.mydata.Variables.(char(plotvar)).units, 'fontsize', 16);
    end


    if cstlns == 1
        if ~isfield(handles, 'cstlns')
            load coast
            handles.cstlns = [long lat];
            % Update handles
            guidata(hObject, handles);
        end
        plot(handles.cstlns(:, 1), handles.cstlns(:, 2), 'k', 'linewidth', 2)
    end
    
    if ctchmnts == 1
        if ~isfield(handles, 'ctchlns')
            load ctchmnt_bnds_grdc.mat
            handles.ctchlns = ctchmnt_bnds;
            % Update handles
            guidata(hObject, handles);
        end
        plot(handles.ctchlns(:, 1), handles.ctchlns(:, 2), 'b', 'linewidth', 2)
    end
    
    if brdrs == 1
        if ~isfield(handles, 'brdrs')
            load natbrdrs_TMWB3.mat
            handles.brdrs = brdrs;
            % Update handles
            guidata(hObject, handles);
        end
        plot(handles.brdrs(:, 1), handles.brdrs(:, 2), 'r', 'linewidth', 2)
    end
 
    % Flip Y-Axis
    axis xy

    % Further options
    % Set the selected colormap
    colormap(clrmps) 
    % Set the selected caxis
    caxis(cxis)      
    % Truncate the spatial extent map according to the data
    axis([min(lons(:)) max(lons(:)) min(lats(:)) max(lats(:))]) 
  
    % Hold off for new maps!
    hold off
    
else
    warning('off', 'backtrace')
    warning('plotstructmap.m: Reached end of timeseries!')
end



% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
