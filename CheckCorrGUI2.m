function varargout = CheckCorrGUI2(varargin)
% CHECKCORRGUI2 MATLAB code for CheckCorrGUI2.fig
%      CHECKCORRGUI2, by itself, creates a new CHECKCORRGUI2 or raises the existing
%      singleton*.
%
%      H = CHECKCORRGUI2 returns the handle to a new CHECKCORRGUI2 or the handle to
%      the existing singleton*.
%
%      CHECKCORRGUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKCORRGUI2.M with the given input arguments.
%
%      CHECKCORRGUI2('Property','Value',...) creates a new CHECKCORRGUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheckCorrGUI2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheckCorrGUI2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheckCorrGUI2

% Last Modified by GUIDE v2.5 29-Nov-2016 15:23:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheckCorrGUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @CheckCorrGUI2_OutputFcn, ...
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


% --- Executes just before CheckCorrGUI2 is made visible.
function CheckCorrGUI2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheckCorrGUI2 (see VARARGIN)

% Choose default command line output for CheckCorrGUI2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% grab guiData from varargin
guiData = varargin{1};
phantomData = guiData{1}; % phantom image data
searchRegionData = guiData{2}; % search region image data
markerDataPass = guiData{3}; % markers of sphere locations for fusing to the image
markerDataFail = guiData{4}; % markers of sphere locations for fusing to the image
grndTruth = guiData{5}; % ground truth data
% store the data so that it can be accessed by the other functions
setappdata(hObject, 'phantomData', phantomData);
setappdata(hObject, 'searchRegionData', searchRegionData);
setappdata(hObject, 'markerDataPass', markerDataPass);
setappdata(hObject, 'markerDataFail', markerDataFail);
setappdata(hObject, 'grndTruth', grndTruth);
numOfImg = length(phantomData(1,1,:));
setappdata(hObject, 'numOfImg', numOfImg);
% grab the handle of the first plot
axes(handles.Plot1)
% plot the middle slice
sphereDim = length(searchRegionData(1,1,:));
sphereCenterSlice = round(median(1:1:sphereDim)); % find the index of the center sphere slice
setappdata(gcf, 'slice',sphereCenterSlice);
% create green markers for location of the spheres
phantomI = phantomData(:,:,sphereCenterSlice); % the phantom image
markerIPass = markerDataPass(:,:,sphereCenterSlice); % the marker image of sphere fit
markerIFail = markerDataFail(:,:,sphereCenterSlice); % the marker image of sphere fit
grndTruthI = grndTruth(:,:,sphereCenterSlice); % the marker image grnd truth
searchRegionI = searchRegionData(:,:,sphereCenterSlice); 
% fit marker plotting
imshow(phantomI,[])

rgb= [0,183,229]/255;
blue = cat(3, rgb(1).*ones(size(searchRegionI)),rgb(2).*ones(size(searchRegionI)),rgb(3).*ones(size(searchRegionI)));
hold on
hBlue = imshow(blue);
hold off
set(hBlue, 'AlphaData', searchRegionI) % make color sheet only show markers

% cat(3,r,g,b)
green = cat(3, zeros(size(markerIPass)),ones(size(markerIPass)), zeros(size(markerIPass)));
hold on 
hGreen = imshow(green); 
hold off
set(hGreen, 'AlphaData', markerIPass) % make color sheet only show markers

red = cat(3, ones(size(markerIFail)),zeros(size(markerIFail)), zeros(size(markerIFail)));
hold on 
hRed = imshow(red); 
hold off
set(hRed, 'AlphaData', markerIFail) % make color sheet only show markers

% ground truth marker plotting
% cat(3,r,g,b)
white = cat(3, ones(size(grndTruthI)),ones(size(grndTruthI)), ones(size(grndTruthI)));
hold on 
h = imshow(white); 
hold off
set(h, 'AlphaData', grndTruthI) % make color sheet only show markers


title('Spatial Integrity Phantom','FontSize',14)


% UIWAIT makes CheckCorrGUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function updatePlot(handles)
phantomData = getappdata(gcf, 'phantomData');
searchRegionData = getappdata(gcf, 'searchRegionData');
markerDataPass = getappdata(gcf, 'markerDataPass');
markerDataFail = getappdata(gcf, 'markerDataFail');
grndTruth = getappdata(gcf, 'grndTruth');
slice = getappdata(gcf, 'slice');
if slice < 1
    slice = 1; % there is no zero slice
end
% grab the handle of the first plot
axes(handles.Plot1)
currentAxes1 = gca;
currentXLim1 = currentAxes1.XLim;
currentYLim1 = currentAxes1.YLim;
% plot the current slice
% create green markers for location of the spheres
phantomI = phantomData(:,:,slice); % the phantom image
markerIPass = markerDataPass(:,:,slice); % the marker image of sphere fit
markerIFail = markerDataFail(:,:,slice); % the marker image of sphere fit
grndTruthI = grndTruth(:,:,slice); % the marker image grnd truth
searchRegionI = searchRegionData(:,:,slice); 

% fit marker plotting
imshow(phantomI,[])

rgb= [0,183,229]/255;
blue = cat(3, rgb(1).*ones(size(searchRegionI)),rgb(2).*ones(size(searchRegionI)),rgb(3).*ones(size(searchRegionI)));
hold on
hBlue = imshow(blue);
hold off
set(hBlue, 'AlphaData', searchRegionI) % make color sheet only show markers

% cat(3,r,g,b)
green = cat(3, zeros(size(markerIPass)),ones(size(markerIPass)), zeros(size(markerIPass)));
hold on 
hGreen = imshow(green); 
hold off
set(hGreen, 'AlphaData', markerIPass) % make color sheet only show markers

red = cat(3, ones(size(markerIFail)),zeros(size(markerIFail)), zeros(size(markerIFail)));
hold on 
hRed = imshow(red); 
hold off
set(hRed, 'AlphaData', markerIFail) % make color sheet only show markers

% ground truth marker plotting
% cat(3,r,g,b)
white = cat(3, ones(size(grndTruthI)),ones(size(grndTruthI)), ones(size(grndTruthI)));
hold on 
h = imshow(white); 
hold off
set(h, 'AlphaData', grndTruthI) % make color sheet only show markers
title('Spatial Integrity Phantom','FontSize',14)
xlim(currentXLim1)
ylim(currentYLim1)
title('Sphere Model','FontSize',14)

% --- Outputs from this function are returned to the command line.
function varargout = CheckCorrGUI2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
numOfImg =  getappdata(gcf, 'numOfImg');
set(hObject, 'SliderStep', [1/numOfImg , 10/numOfImg ]);
slice=round(get(handles.slider1, 'Value')*numOfImg);

if slice < 1
    slice = 1; % there is no zero slice
end
if slice > numOfImg
   slice = numOfImg; 
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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
searchRegionData = getappdata(gcf, 'searchRegionData');
numOfImg =  getappdata(gcf, 'numOfImg');
slice=round(str2num(get(handles.edit1, 'String')));


if slice > numOfImg
    % value specified larger than number of slices acquired
    slice=round(get(handles.slider1, 'Value')*numOfImg);
    set(handles.edit1, 'String', slice);
else
    % value within bounds of slices, change value of slider
    set(handles.slider1, 'Value', slice/numOfImg);
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


phantomData = getappdata(gcf, 'phantomData');
searchRegionData = getappdata(gcf, 'searchRegionData');
markerDataPass = getappdata(gcf, 'markerDataPass');
markerDataFail = getappdata(gcf, 'markerDataFail');
grndTruth = getappdata(gcf, 'grndTruth');
slice = getappdata(gcf, 'slice');
if slice < 1
    slice = 1; % there is no zero slice
end
% grab the handle of the first plot
axes(handles.Plot1)
% plot the current slice
% create green markers for location of the spheres
phantomI = phantomData(:,:,slice); % the phantom image
markerIPass = markerDataPass(:,:,slice); % the marker image of sphere fit
markerIFail = markerDataFail(:,:,slice); % the marker image of sphere fit
grndTruthI = grndTruth(:,:,slice); % the marker image grnd truth
searchRegionI = searchRegionData(:,:,slice); 
% fit marker plotting
imshow(phantomI,[])

rgb= [0,183,229]/255;
blue = cat(3, rgb(1).*ones(size(searchRegionI)),rgb(2).*ones(size(searchRegionI)),rgb(3).*ones(size(searchRegionI)));
hold on
hBlue = imshow(blue);
hold off
set(hBlue, 'AlphaData', searchRegionI) % make color sheet only show markers

% cat(3,r,g,b)
green = cat(3, zeros(size(markerIPass)),ones(size(markerIPass)), zeros(size(markerIPass)));
hold on 
hGreen = imshow(green); 
hold off
set(hGreen, 'AlphaData', markerIPass) % make color sheet only show markers

red = cat(3, ones(size(markerIFail)),zeros(size(markerIFail)), zeros(size(markerIFail)));
hold on 
hRed = imshow(red); 
hold off
set(hRed, 'AlphaData', markerIFail) % make color sheet only show markers

% ground truth marker plotting
% cat(3,r,g,b)
white = cat(3, ones(size(grndTruthI)),ones(size(grndTruthI)), ones(size(grndTruthI)));
hold on 
h = imshow(white); 
hold off
set(h, 'AlphaData', grndTruthI) % make color sheet only show markers
title('Spatial Integrity Phantom','FontSize',14)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
