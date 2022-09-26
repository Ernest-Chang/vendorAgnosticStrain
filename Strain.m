function varargout = Strain(varargin)
% STRAIN MATLAB code for Strain.fig
%      STRAIN, by itself, creates a new STRAIN or raises the existing
%      singleton*.
%
%      H = STRAIN returns the handle to a new STRAIN or the handle to
%      the existing singleton*.
%
%      STRAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRAIN.M with the given input arguments.
%
%      STRAIN('Property','Value',...) creates a new STRAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Strain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Strain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Strain

% Last Modified by GUIDE v2.5 20-Nov-2021 16:28:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Strain_OpeningFcn, ...
                   'gui_OutputFcn',  @Strain_OutputFcn, ...
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


% --- Executes just before Strain is made visible.
function Strain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Strain (see VARARGIN)

% Choose default command line output for Strain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% h=findobj('Tag','pushbutton1');
% dataStruc=get(h,'UserData');
% ampMap=dataStruc.ampMap;
% aMap=dataStruc.aMap;
% nF=dataStruc.nF;


% UIWAIT makes Strain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Strain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hData, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 
[fileName,pathName]=uigetfile('*.*');
info =VideoReader([pathName fileName]);
nF=round(info.NumFrames);
Y1=read(info);
Y2=squeeze(Y1(:,:,1,:)/3+Y1(:,:,2,:)/3+Y1(:,:,3,:)/3)*1.3;

axes(handles.axes3);
imagesc(Y2(:,:,1)),colormap(gray);

data=struct('Y',Y2,'nF',nF,'pname',pathName,'fname',fileName);
set(hData,'UserData',data);

h=findobj('Tag','slider1');
set(h,'Value',round(nF/2));
set(h, 'Max', nF);
set(h, 'Min', 1);

h=findobj('Tag','slider2');
set(h,'Value',round(nF/2));
set(h, 'Max', nF);
set(h, 'Min', 1);



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(haxes1, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(haxes2, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3

% --- Executes on slider movement.
function slider1_Callback(hSlider, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%retrieve data from pushbutton1
h=findobj('Tag','pushbutton1');
dataStruc=get(h,'UserData');

rawData=dataStruc.Y;
noF=dataStruc.nF;

% 
% set(hSlider,'Value',1);
set(hSlider, 'Min', 1);
set(hSlider, 'Max', noF);


sliderVal=round(get(hSlider,'Value'));
data=struct('Y',rawData,'nF',noF,'init',sliderVal);
set(hSlider,'UserData',data);

display(sliderVal)

axes(handles.axes3);
imagesc(rawData(:,:,sliderVal)),colormap(gray)

h=findobj('Tag','slider2');
set(h,'Value',sliderVal);
% 
% hline1 = imdistline(gca);
% hline2 = imdistline(gca);
% hline3 = imdistline(gca);
% set(h, 'Max', nF);
% set(h, 'Min', 1);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hSlider1, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hSlider1,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hSlider1,'BackgroundColor',[.9 .9 .9]);
end
% 







% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hExecute, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


mymap=[
             0         0         0
    0.1250         0         0
    0.2500         0         0
    0.3750         0         0
    0.5000         0         0
    0.6250         0         0
    0.7500         0         0
    0.8750         0         0
    1.0000         0         0
    1.0000    0.0625         0
    1.0000    0.1250         0
    1.0000    0.1875         0
    1.0000    0.2500         0
    1.0000    0.3125         0
    1.0000    0.3750         0
    1.0000    0.4375         0
    1.0000    0.5000         0
    1.0000    0.5625         0
    1.0000    0.6250         0
    1.0000    0.6875         0
    1.0000    0.7500         0
    1.0000    0.8125         0
    1.0000    0.8750         0
    1.0000    0.9375         0
    1.0000    1.0000         0
    0.9375    1.0000    0.0625
    0.8750    1.0000    0.1250
    0.8125    1.0000    0.1875
    0.7500    1.0000    0.2500
    0.6875    1.0000    0.3125
    0.6250    1.0000    0.3750
    0.5625    1.0000    0.4375
    0.5000    1.0000    0.5000
    0.4375    1.0000    0.5625
    0.3750    1.0000    0.6250
    0.3125    1.0000    0.6875
    0.2500    1.0000    0.7500
    0.1875    1.0000    0.8125
    0.1250    1.0000    0.8750
    0.0625    1.0000    0.9375
         0    1.0000    1.0000
         0    0.9375    1.0000
         0    0.8750    1.0000
         0    0.8125    1.0000
         0    0.7500    1.0000
         0    0.6875    1.0000
         0    0.6250    1.0000
         0    0.5625    1.0000
         0    0.5000    1.0000
         0    0.4375    1.0000
         0    0.3750    1.0000
         0    0.3125    1.0000
         0    0.2500    1.0000
         0    0.1875    1.0000
         0    0.1250    1.0000
         0    0.0625    1.0000
         0         0    1.0000
         0         0    0.9375
         0         0    0.8750
         0         0    0.8125
         0         0    0.7500
         0         0    0.6875
         0         0    0.6250
         0         0    0.5625];

     
     
h=findobj('Tag','slider1');
dataStruc=get(h,'UserData');
init=dataStruc.init;
nF=dataStruc.nF;
data=dataStruc.Y;

     
h=findobj('Tag','slider2');
dataStruc=get(h,'UserData');
fin=dataStruc.fin;

h=findobj('Tag','pushbutton1');
dataStruc=get(h,'UserData');
path=dataStruc.pname;
file=dataStruc.fname;


[ampMap,aMap,nF,mStrain,cycleL]=speckle_strain(init,fin,nF,data,path,file);

     
axes(handles.axes1);
imagesc(ampMap),colormap(jet),caxis([0 0.40]), title('mean strain is ',-mStrain*100),colorbar
% F = getframe(handles.axes1);
% Image = frame2im(F);

newfig=figure;
imagesc(ampMap),colormap(jet),caxis([0 0.40]), title('mean strain is ',-mStrain*100),colorbar
newimg=frame2im(getframe(newfig));
imwrite(newimg, [path file(1:end-4) '.bmp'])



axes(handles.axes2);
imagesc(aMap),colormap(mymap), caxis([0 cycleL]),colormap(mymap),title('Activation map'),colorbar;
F = getframe(handles.axes2);
Image2 = frame2im(F);
imwrite(Image2, 'aMap.jpg')

newfig2=figure;
imagesc(medfilt2(aMap,[3 3])),colormap(mymap),caxis([0 cycleL]),colorbar
newimg=frame2im(getframe(newfig2));
imwrite(newimg, [path file(1:end-4) 'TimeMap.bmp'])


% --- Executes on slider movement.
function slider2_Callback(hSlider2, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%retrieve data from pushbutton1
h=findobj('Tag','pushbutton1');
dataStruc=get(h,'UserData');

rawData=dataStruc.Y;
noF=dataStruc.nF;

% 
% set(hSlider,'Value',1);
set(hSlider2, 'Min', 1);
set(hSlider2, 'Max', noF);


sliderVal=round(get(hSlider2,'Value'));
data=struct('fin',sliderVal);
set(hSlider2,'UserData',data);

% display(sliderVal)

axes(handles.axes3);
imagesc(rawData(:,:,sliderVal)),colormap(gray)





% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
