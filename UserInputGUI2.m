function varargout = UserInputGUI2(varargin)
% USERINPUTGUI2 MATLAB code for UserInputGUI2.fig
%      USERINPUTGUI2, by itself, creates a new USERINPUTGUI2 or raises the existing
%      singleton*.
%
%      H = USERINPUTGUI2 returns the handle to a new USERINPUTGUI2 or the handle to
%      the existing singleton*.
%
%      USERINPUTGUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in USERINPUTGUI2.M with the given input arguments.
%
%      USERINPUTGUI2('Property','Value',...) creates a new USERINPUTGUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UserInputGUI2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UserInputGUI2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UserInputGUI2

% Last Modified by GUIDE v2.5 29-Sep-2016 15:54:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UserInputGUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @UserInputGUI2_OutputFcn, ...
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


% --- Executes just before UserInputGUI2 is made visible.
function UserInputGUI2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UserInputGUI2 (see VARARGIN)

% Choose default command line output for UserInputGUI2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% grab guiData from varargin
guiData = varargin{1};
fileData = guiData{1}; % image data
fileMap = guiData{2}; % image map

setappdata(gcf, 'fileData', fileData);
setappdata(gcf, 'fileMap', fileMap);


% grab the handle of the first plot
axes(handles.Plot1)

% plot the middle slice
midSlice = round(length(fileData(1,1,:))/2);
imshow(fileData(:,:,midSlice),fileMap{midSlice})
title('Spatial Integrity Phantom','FontSize',14)

setappdata(gcf, 'slice',midSlice);
% UIWAIT makes UserInputGUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function updatePlot(handles)
fileData = getappdata(gcf, 'fileData');
fileMap = getappdata(gcf, 'fileMap');
slice = getappdata(gcf, 'slice');

% grab the handle of the first plot
axes(handles.Plot1)
currentAxes = gca;
currentXLim = currentAxes.XLim;
currentYLim = currentAxes.YLim;
% plot the current slice
imshow(fileData(:,:,slice),fileMap{slice})
xlim(currentXLim)
ylim(currentYLim)

title('Spatial Integrity Phantom','FontSize',14)

% --- Outputs from this function are returned to the command line.
function varargout = UserInputGUI2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slice=round(str2num(get(handles.edit1, 'String')));
fileData = getappdata(gcf, 'fileData');

if slice > length(fileData(1,1,:))
    % value specified larger than number of slices acquired
    slice=round(get(handles.slider1, 'Value')*length(fileData(1,1,:)));
    set(handles.edit1, 'String', slice);
else
    % value within bounds of slices, change value of slider
    set(handles.slider1, 'Value', slice/length(fileData(1,1,:)));
    setappdata(gcf, 'slice', slice);
end
updatePlot(handles);
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


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileData = getappdata(gcf, 'fileData');
numOfImg = length(fileData(1,1,:));
set(hObject, 'SliderStep', [1/numOfImg , 10/numOfImg ]);
slice=round(get(handles.slider1, 'Value')*length(fileData(1,1,:)));
if slice < 1
    slice = 1; % there is no zero slice
end
setappdata(gcf, 'slice', slice);
set(handles.edit1, 'String', num2str(slice));
updatePlot(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
fileData = getappdata(gcf, 'fileData');
fileMap = getappdata(gcf, 'fileMap');
slice = getappdata(gcf, 'slice');

% grab the handle of the first plot
axes(handles.Plot1);

% plot the current slice
imshow(fileData(:,:,slice),fileMap{slice})
title('Spatial Integrity Phantom','FontSize',14)