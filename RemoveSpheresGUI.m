function varargout = RemoveSpheresGUI(varargin)
% REMOVESPHERESGUI MATLAB code for RemoveSpheresGUI.fig
%      REMOVESPHERESGUI, by itself, creates a new REMOVESPHERESGUI or raises the existing
%      singleton*.
%
%      H = REMOVESPHERESGUI returns the handle to a new REMOVESPHERESGUI or the handle to
%      the existing singleton*.
%
%      REMOVESPHERESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REMOVESPHERESGUI.M with the given input arguments.
%
%      REMOVESPHERESGUI('Property','Value',...) creates a new REMOVESPHERESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RemoveSpheresGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RemoveSpheresGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RemoveSpheresGUI

% Last Modified by GUIDE v2.5 30-Nov-2016 11:21:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RemoveSpheresGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RemoveSpheresGUI_OutputFcn, ...
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


% --- Executes just before RemoveSpheresGUI is made visible.
function RemoveSpheresGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RemoveSpheresGUI (see VARARGIN)

% Choose default command line output for RemoveSpheresGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% clear any stored sphere locations that need to be removed from previous
% iterations of calling this GUI
if isappdata(hObject,'spheresToRemove')
    rmappdata(hObject,'spheresToRemove')
end
% grab guiData from varargin
guiData = varargin{1};
phantomData = guiData{1}; % phantom image data
searchRegionData = guiData{2}; % search region image data
markerDataPass = guiData{3}; % markers of sphere locations for fusing to the image
markerDataFail = guiData{4}; % markers of sphere locations for fusing to the image
grndTruth = guiData{5}; % ground truth data

% create markers that are going to be removed
volRemove = zeros(length(phantomData(:,1,1)),length(phantomData(1,:,1)),...
    length(phantomData(1,1,:)));

% store the data so that it can be accessed by the other functions
setappdata(hObject, 'countSpheresToRemove',0)
setappdata(hObject, 'phantomData', phantomData);
setappdata(hObject, 'searchRegionData', searchRegionData);
setappdata(hObject, 'markerDataPass', markerDataPass);
setappdata(hObject, 'markerDataFail', markerDataFail);
setappdata(hObject, 'grndTruth', grndTruth);
setappdata(hObject, 'volRemove',volRemove);
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
removeI = volRemove(:,:,sphereCenterSlice);

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

red =  cat(3, 255.*ones(size(removeI)),0.*ones(size(removeI)), zeros(size(removeI)));
hold on 
h = imshow(red); 
hold off
set(h, 'AlphaData', removeI) % make color sheet only show markers

title('Spatial Integrity Phantom','FontSize',14)


% UIWAIT makes RemoveSpheresGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);

function updatePlot(handles)
phantomData = getappdata(gcf, 'phantomData');
searchRegionData = getappdata(gcf, 'searchRegionData');
markerDataPass = getappdata(gcf, 'markerDataPass');
markerDataFail = getappdata(gcf, 'markerDataFail');
grndTruth = getappdata(gcf, 'grndTruth');
slice = getappdata(gcf, 'slice');
volRemove = getappdata(gcf,'volRemove');
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
removeI = volRemove(:,:,slice);

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

red =  cat(3, 255.*ones(size(removeI)),0.*ones(size(removeI)), zeros(size(removeI)));
hold on 
h = imshow(red); 
hold off
set(h, 'AlphaData', removeI) % make color sheet only show markers

title('Spatial Integrity Phantom','FontSize',14)
xlim(currentXLim1)
ylim(currentYLim1)
title('Sphere Model','FontSize',14)

% --- Outputs from this function are returned to the command line.
function varargout = RemoveSpheresGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if isequal(get(hObject, 'waitstatus'), 'waiting')
% countSpheresToRemove = getappdata(hObject,'countSpheresToRemove');
% else
%     
% end
% if countSpheresToRemove > 0
%     varargout{1} = getappdata(hObject,'spheresToRemove');
% else
%     varargout{1} = NaN;
% end
% uiresume(handles.figure1);
% Get default command line output from handles structure
varargout{1} = getappdata(gcf,'spheresToRemove');

delete(handles.figure1)


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
volRemove = getappdata(gcf,'volRemove');

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
removeI = volRemove(:,:,slice);
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


red =  cat(3, 255.*ones(size(removeI)),0.*ones(size(removeI)), zeros(size(removeI)));
hold on 
h = imshow(red); 
hold off
set(h, 'AlphaData', removeI) % make color sheet only show markers


title('Spatial Integrity Phantom','FontSize',14)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
countSpheresToRemove = getappdata(gcf,'countSpheresToRemove');
volRemove = getappdata(gcf,'volRemove');
% if this isn't the first sphere to be removed obtain the others in the
% list
if countSpheresToRemove ~= 0
    spheresToRemove = getappdata(gcf,'spheresToRemove');
end
countSpheresToRemove = countSpheresToRemove + 1;
setappdata(gcf,'countSpheresToRemove',countSpheresToRemove)
% obtain the data locations
[xRemove, yRemove] = ginput(1);
zRemove = getappdata(gcf,'slice');
spheresToRemove(countSpheresToRemove,1) = xRemove;
spheresToRemove(countSpheresToRemove,2) = yRemove;
spheresToRemove(countSpheresToRemove,3) = zRemove;
% calculate a little square for marking the spheres to be removed
xDim = length(volRemove(1,:,1));
yDim = length(volRemove(:,1,1));
zDim = length(volRemove(1,1,:));
xPos = round(xRemove);
yPos = round(yRemove);
zPos = round(zRemove);
% ensure location does not exceed boundary
% y-location
if (yPos > 1)&&(yPos < yDim)
    yLoc = (yPos-1):1:(yPos+1);
end
if (yPos <= 1)
    yLoc = 1:1:(2);
end
if (yPos >= yDim)
    yLoc = (yDim-1):1:yDim;
end
% x-location
if (xPos > 1)&&(xPos < xDim)
    xLoc = (xPos-1):1:(xPos+1);
end
if (xPos <= 1)
    xLoc = 1:1:(2);
end
if (xPos >= xDim)
    xLoc = (xDim-1):1:xDim;
end
% z-location
if (zPos > 1)&&(zPos < zDim)
    zLoc = (zPos-1):1:(zPos+1);
end
if (zPos <= 1)
    zLoc = 1:1:(2);
end
if (zPos >= zDim)
    zLoc = (zDim-1):1:zDim;
end
volRemove(yLoc,xLoc,zLoc) = 1;
setappdata(gcf,'spheresToRemove',spheresToRemove);
setappdata(gcf,'volRemove',volRemove)
updatePlot(handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
