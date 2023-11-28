function varargout = sem_gui(varargin)
% SEM_GUI M-file for sem_gui.fig
%      SEM_GUI, by itself, creates a new SEM_GUI or raises the existing
%      singleton*.
%
%      H = SEM_GUI returns the handle to a new SEM_GUI or the handle to
%      the existing singleton*.
%
%      SEM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEM_GUI.M with the given input arguments.
%
%      SEM_GUI('Property','Value',...) creates a new SEM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sem_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sem_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sem_gui

% Last Modified by GUIDE v2.5 29-Nov-2012 12:58:33

% Begin initialization code - DO NOT EDIT

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sem_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sem_gui_OutputFcn, ...
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


% --- Executes just before sem_gui is made visible.
function sem_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sem_gui (see VARARGIN)

%clear all;

% initialising
model_specs=struct([]);
data=[];
pars=[];
minval=[];
gof_measures=struct([]);
std_res=[];
se=[];
ci=[];
se_bt=[];
ci_bt=[];

assignin('base','model_specs',model_specs);
assignin('base','data',data);
assignin('base','pars',pars);
assignin('base','minval',minval);
assignin('base','gof_measures',gof_measures);
assignin('base','std_res',std_res);
assignin('base','se',se);
assignin('base','ci',ci);
assignin('base','se_bt',se_bt);
assignin('base','ci_bt',ci_bt);

% Choose default command line output for sem_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sem_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sem_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
num_obs = str2double(get(hObject,'string'));
if isnan(num_obs)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
  return;
end
% Proceed with callback...


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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
num_lat = str2double(get(hObject,'string'));
if isnan(num_lat)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
	return
end
% Proceed with callback...


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
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

% get the dimension of the symmetric and asymmetric matrices

num_obs=str2double(get(handles.edit2,'String'));
num_lat=str2double(get(handles.edit3,'String'));
var_num=num_obs+num_lat;

sym_data=get(handles.sym_matrix,'Data');
asym_data=get(handles.asym_matrix,'Data');

% checking for non-reciprocal symmetric links

for i=1:var_num,
    
    for j=1:var_num,
        
        if i~=j
            
           if sym_data(i,j)~=sym_data(j,i)
               
              errordlg('Two-way paths should be reflected on both sides of the diagonal','Bad Input','modal');
              uicontrol(hObject); 
              return; 
              
           end
            
            
        end
        
    end
    
end

% checking whether variances of observed variables are free parameters

for i=1:var_num
   
     if ((sym_data(i,i)>0)==0),  
     errordlg('Variances of observed variables should be >0 or free parameters','Bad Input','modal');
     uicontrol(hObject); 
     return;    
        
    end
    
end

% checking for reciprocal asymmetric links

for i=1:var_num,
    
    for j=1:var_num,
        
        if i~=j
            
           if (asym_data(i,j)>0) && (asym_data(j,i)>0)
               
              errordlg('One-way paths cannot be replicated on both sides of the diagonal','Bad Input','modal');
              uicontrol(hObject); 
              return; 
              
           end
            
            
        end
        
    end
    
end

model_specs=struct('asym_data',asym_data,'sym_data',sym_data,'obs',num_obs,'lat',num_lat);

assignin('base','model_specs',model_specs);

% disabling text boxes and Create button

set(handles.pushbutton8,'Enable','Off');

% re-sizing and re-titling panels

set(handles.uipanel2,'Units','normalized','Position',[60/857 230/689 400/857 330/689],'Title','Parameter details','FontUnits','normalized');
set(handles.uipanel3,'Units','normalized','Position',[60/857 80/689 400/857 120/689],'Title','Fit details','FontUnits','normalized');

% activating model estimation buttons

set(handles.pushbutton2,'Enable','On');
set(handles.pushbutton9,'Enable','On');
set(handles.pushbutton10,'Enable','On');
set(handles.popupmenu1,'Enable','On');

% deleting components and invisibling components

delete(handles.sym_matrix);
delete(handles.asym_matrix);
set(hObject,'Visible','Off');

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
error_flag_pars=0;
error_flag_id=0;

% imports model details from base workspace
model=evalin('base','model_specs');

% calculates number of free parameters
symel=0;
for idx1=1:size(model.sym_data,1),
    
    for idx2=1:idx1,
        
        if model.sym_data(idx1,idx2)==Inf,
   
        symel=symel+1;
        
        end
        
    end
    
end

par_num=symel+length(find(model.asym_data==Inf))-model.lat;

% Input dialog box for initial parameter values
prompt = {['Enter initial guess - ',int2str(par_num),' values:']};
dlg_title = 'Initial parameter values';
num_lines = par_num;
def={[char(num2str(ones(par_num,1)))]};
answer = inputdlg(prompt,dlg_title,num_lines,def);
val = str2num(answer{1}); %#ok<ST2NM>

[l,b]=size(val);

if l~=par_num,
   
  errordlg('Incorrect number of parameter values - please re-enter values','Initial parameter values','modal');
  error_flag_pars=1;
    
end

if error_flag_pars==0

% checking for possibility of identification
num_nr=(model.obs*(model.obs+1))/2;
if num_nr<par_num
    
    errordlg('Model cannot be identified - number of unknowns (parameters) exceeds number of knowns (non-redundant elements)','Invalid model','modal');
    error_flag_id=1;
    
end

if error_flag_id==0
    
% transposing parameter values

val=val';

% estimating parameters
option = get(handles.popupmenu1,'Value');

switch option
    
case 1  
 [pars,minval]=anneal(@model_est_ml,val);   
case 2
 [pars,minval]=anneal(@model_est_gls,val);
case 3      
 [pars,minval]=anneal(@model_est_uls,val);
end
 
% writing output to base workspace
assignin('base','pars',pars);
assignin('base','minval',minval);

% Tabulate data
columns={'Parameters','SE','CI-low 95%','CI-high 95%','SE(BT)','CI(BT)-low 95%','CI(BT)-high 95%'};
table_data=horzcat(pars');

position=[80/857 240/689 140/857 300/689];
t=uitable('Data',table_data,'ColumnName',columns,'Parent',handles.figure1,'Units','normalized','Position',position,'ColumnWidth','auto','FontUnits','normalized','FontSize',0.043);

% adding handle to table
handles.t=t;
guidata(hObject,handles);

% checking ratio of samples to free parameters
data=evalin('base','data');
samp_ratio=size(data,1)/par_num;

if samp_ratio<10   
    warndlg('Number of samples is likely too few for given number of parameters','Low sample size','modal');
end

end

end

% enabling relevant buttons

set(handles.pushbutton3,'Enable','On')
set(handles.pushbutton11,'Enable','On')
set(handles.pushbutton15,'Enable','On')



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% estimating parameters
option = get(handles.popupmenu1,'Value');

    switch option
    
    case 1 
    
    pars=evalin('base','pars');
    minval=evalin('base','minval');

    % calculate confidence intervals
    [ci,se]=model_ci(pars,minval);
    assignin('base','ci',ci);
    assignin('base','se',se);

    % Tabulate data
    columns={'Parameters','SE','CI-low 95%','CI-high 95%','SE(BT)','CI(BT)-low 95%','CI(BT)-high 95%'};
    table_data=horzcat(pars',se',ci(2,:)',ci(1,:)');
    position=[80/857 240/689 360/857 300/689];
    set(handles.t,'Data',table_data,'ColumnName',columns,'Parent',handles.figure1,'Units','normalized','Position',position,'ColumnWidth','auto','FontUnits','normalized','FontSize',0.043);
    
    %Enabling relevant buttons
    set(handles.pushbutton4,'Enable','On');
    set(handles.pushbutton5,'Enable','On');
    set(handles.pushbutton6,'Enable','On');
    set(handles.pushbutton7,'Enable','On');

    case 2
        
        h=msgbox('SEs are computed only for the ML option','SE estimation','modal');
        
    case 3
        
        h=msgbox('SEs are computed only for the ML option','SE estimation','modal');
        
    end
        
        

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% getting parameter values and minimum values
pars=evalin('base','pars');
minval=evalin('base','minval');

% goodness-of-fit measures
option = get(handles.popupmenu1,'Value');

switch option
    
case 1    
gof_measures=goftests_ml(pars,minval); 
case 2
gof_measures=goftests_gls(pars,minval);
case 3       
gof_measures=goftests_uls(pars,minval);
end

% writing goodness-of-fit measures to workspace
assignin('base','gof_measures',gof_measures);

% Tabulate data
columns={'Chi-square','NFI','NNFI','AIC','SABIC','RMSEA'};

for i=1:3,
    
   gof_data(i,1)=gof_measures(1,i).CHISQ;
   gof_data(i,2)=gof_measures(1,i).NFI;
   gof_data(i,3)=gof_measures(1,i).NNFI;
   gof_data(i,4)=gof_measures(1,i).AIC;
   gof_data(i,5)=gof_measures(1,i).SABIC;
   gof_data(i,6)=gof_measures(1,i).RMSEA;
    
end

position=[70/857 90/689 380/857 90/689];
gof_table=uitable('Data',gof_data,'ColumnName',columns,'Units','normalized','Position',position,'ColumnWidth',{58},'FontUnits','normalized','FontSize',0.15);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

error_flag=0;

model=evalin('base','model_specs');
num_obs=model.obs;

% Input dialog box for initial parameter values
prompt = {['Enter ',int2str(num_obs),' labels:']};
dlg_title = 'Variable labels';
num_lines = num_obs;
def={[char(num2str([1:num_obs]'))]};
response = inputdlg(prompt,dlg_title,num_lines,def);

[l,b]=size(response{1,1});

% checking if number of labels is correct

if l~=num_obs,
   
  errordlg('Incorrect number of labels - please re-enter labels','Label details','modal');
  error_flag=1;
    
end

if error_flag==0
    
pars=evalin('base','pars');
[raw_res,std_res]=residuals(pars);

% Plotting image of residuals

figure;
imagesc(std_res);
colorbar;
title('Standardised residuals');
xlabel('Observed variables');
ylabel('Observed variables');
set(gca,'XTick',[1:num_obs]','XTickLabel',response{1});
set(gca,'YTick',[1:num_obs]','YTickLabel',response{1});

% writing to workspace
assignin('base','std_res',std_res);

end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject -  handle to pushbutton8 (see GCBO)
% eventdata - reserved - to be defined in a future version of MATLAB
% handles - structure with handles and user data (see GUIDATA)
num_obs=str2double(get(handles.edit2,'String'));
num_lat=str2double(get(handles.edit3,'String'));

var_num=num_obs+num_lat;

sym_data=zeros(var_num,var_num);
asym_data=zeros(var_num,var_num);

% Tabulate data

position=[15/470 10/250 440/470 230/250];
sym_matrix=uitable('Data',sym_data,'Parent',handles.uipanel2,'Units','normalized','Position',position,'ColumnWidth',{35},'ColumnEditable',logical(ones(1,var_num)),'FontUnits','normalized','FontSize',0.062,'CellEditCallback',{@sym_editcback});
asym_matrix=uitable('Data',asym_data,'Parent',handles.uipanel3,'Units','normalized','Position',position,'ColumnWidth',{35},'ColumnEditable',logical(ones(1,var_num)),'FontUnits','normalized','FontSize',0.062,'CellEditCallback',{@asym_editcback});

handles.sym_matrix=sym_matrix;
handles.asym_matrix=asym_matrix;

guidata(hObject,handles);
guidata(hObject,handles);

% enabling Confirm Model button

set(handles.edit2,'Enable','Off');
set(handles.edit3,'Enable','Off');
set(handles.pushbutton1,'Enable','On');


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=uiimport('-file');
eval(['data=s.' char(fieldnames(s)) ';']);
assignin('base','data',data);
%clear('base',char(fieldnames(s)));



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
alpha=0.05;
data=evalin('base','data');

[p1,p1c,p2]=multskewkurt(data,alpha);

% calculating if small or large sample
model=evalin('base','model_specs');

symel=0;
for idx1=1:size(model.sym_data,1),
    
    for idx2=1:idx1,
        
        if model.sym_data(idx1,idx2)==Inf,
   
        symel=symel+1;
        
        end
        
    end
    
end

par_num=symel+length(find(model.asym_data==Inf))-model.lat;
sample_ratio=(size(data,1)/par_num);

% checking if assumptions are met

if p2 <=alpha 
  warndlg('Assumption of multivariate normality is violated (Kurtosis)','Assumption verification','modal');    
end

if (sample_ratio >= 10)

    if  p1 <=alpha
        warndlg('Assumption of multivariate normality is violated (Skewness)','Assumption verification','modal');
    end

    if p1>alpha && p2>alpha
        h=msgbox('Assumption of multivariate normality is valid','Assumption verification','modal');  
    end

else
  
    if p1c <=alpha
        warndlg('Assumption of multivariate normality is violated (Skewness)','Assumption verification','modal');   
    end
    
    if p1c>alpha && p2>alpha
        h=msgbox('       Assumption of multivariate normality is valid','Assumption verification','modal');  
    end 

end




% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

error_flag_id=0;

pars=evalin('base','pars');

% getting number of parameters
num_pars=length(pars);

% estimating the Hessian matrix
option = get(handles.popupmenu1,'Value');
 
switch option
    
case 1   
 [H]=hessian(@model_est_ml,pars);   
case 2
 [H]=hessian(@model_est_gls,pars);
case 3      
 [H]=hessian(@model_est_uls,pars);
end

% estimating the Information matrix
I=-1*H;

% calculating rank of Information matrix
r=rank(I);

% display appropriate message
if r~=num_pars,
    
    errordlg('Model is not identified - please re-specify model','Model ID','modal');
    error_flag_id=1;
    
else
      
  %  message box
   msgbox('          Model is identified','Model ID','modal');
   
end

% Enabling relevant buttons
set(handles.pushbutton4,'Enable','On');
set(handles.pushbutton5,'Enable','On');
set(handles.pushbutton6,'Enable','On');
set(handles.pushbutton7,'Enable','On');

if error_flag_id==1
    
%  re-drawing

delete(handles.t);
set(handles.uipanel2,'Title','Symmetric','Units','pixels','Position',[45 414 621 310]);
set(handles.uipanel3,'Title','Asymmetric','Units','pixels','Position',[45 86 621 310]);
   
model=evalin('base','model_specs');

symel=0;
for idx1=1:size(model.sym_data,1),
    
    for idx2=1:idx1,
        
        if model.sym_data(idx1,idx2)==Inf,
   
        symel=symel+1;
        
        end
        
    end
    
end

par_num=symel+length(find(model.asym_data==Inf))-model.lat;

x=ones(1,par_num);

% calculate F matrix

F=zeros(model.obs,model.obs+model.lat);

for i=1:model.obs,
       
    F(i,model.lat+i)=1;
    
end

% generating symmetric matrix

S=zeros(size(model.sym_data));
S=model.sym_data;
[sidxs(:,1),sidxs(:,2)]=find(model.sym_data==Inf);
sidxsind=find((sidxs(:,2)<=sidxs(:,1))==1);
sidxs=sidxs(sidxsind,:);
snum=length(sidxs);
sidxs=sub2ind(size(S),sidxs(:,1),sidxs(:,2));
S(sidxs)=x(1,1:snum);
for idx1=1:size(S,1),
    for idx2=1:idx1,
        if idx1~=idx2,
        S(idx2,idx1)=S(idx1,idx2);
        end
    end
end

% generating asymmetric matrix

A=zeros(size(model.asym_data));
A=model.asym_data;

for cidx=1:model.lat,

    tek1=find(A(:,cidx)==Inf);
    tek2= find(tek1>model.lat);
    tek1=tek1(tek2,1);
    A(tek1(1,1),cidx)=1;

end

[aidxs]=find(A==Inf);
A(aidxs)=x(1,snum+1:length(x));

position=[40 40 540 230];
sym_matrix=uitable('Data',S,'Parent',handles.uipanel2,'Position',position,'ColumnWidth',{35},'FontSize',10,'CellEditCallback',{@sym_editcback});
asym_matrix=uitable('Data',A,'Parent',handles.uipanel3,'Position',position,'ColumnWidth',{35},'FontSize',10,'CellEditCallback',{@asym_editcback});

handles.sym_matrix=sym_matrix;
handles.asym_matrix=asym_matrix;

guidata(hObject,handles);
guidata(hObject,handles);

set(handles.pushbutton1,'Visible','On');

% de-activating model estimation buttons

set(handles.pushbutton2,'Enable','Off');
set(handles.pushbutton3,'Enable','Off');
set(handles.pushbutton9,'Enable','Off');
set(handles.pushbutton10,'Enable','Off');
set(handles.pushbutton11,'Enable','Off');
set(handles.pushbutton15,'Enable','Off');
set(handles.popupmenu1,'Enable','Off');
    
end

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


function sym_editcback(hObject, eventdata)
    % hObject    Handle to uitable1 (see GCBO)
    % eventdata  Currently selected table indices
    % Callback to erase and replot markers, showing only those
    
 if (eventdata.NewData>=0 && eventdata.NewData<=1) || (eventdata.NewData==Inf)    
     
 else
       tableData = get(hObject, 'data');
       tableData(eventdata.Indices(1), eventdata.Indices(2)) = eventdata.PreviousData;
       set(hObject, 'data', tableData);
       error('Fixed parameters should be between 0 and 1, Free parameters should be Inf');
 end

function asym_editcback(hObject, eventdata)
    % hObject    Handle to uitable1 (see GCBO)
    % eventdata  Currently selected table indices
    % Callback to erase and replot markers, showing only those

% --- Executes when selected cell(s) is changed in uitable6.

if (eventdata.NewData>=0 && eventdata.NewData<=1) || (eventdata.NewData==Inf)    
     
else
       tableData = get(hObject, 'data');
       tableData(eventdata.Indices(1), eventdata.Indices(2)) = eventdata.PreviousData;
       set(hObject, 'data', tableData);
       error('Fixed parameters should be between 0 and 1, Free parameters should be Inf');
 end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

closeGUI = handles.figure1; %handles.figure1 is the GUI figure

guiPosition = get(handles.figure1,'Position'); %get the position of the GUI
guiName = get(handles.figure1,'Name'); %get the name of the GUI
eval(guiName) %call the GUI again

close(closeGUI); %close the old GUI
set(gcf,'Position',guiPosition); %set the position for the new GUI

% re-initialising
model_specs=struct([]);
data=[];
pars=[];
minval=[];
gof_measures=struct([]);
std_res=[];
se=[];
ci=[];
se_bt=[];
ci_bt=[];

assignin('base','model_specs',model_specs);
assignin('base','data',data);
assignin('base','pars',pars);
assignin('base','minval',minval);
assignin('base','gof_measures',gof_measures);
assignin('base','std_res',std_res);
assignin('base','se',se);
assignin('base','ci',ci);
assignin('base','se_bt',se_bt);
assignin('base','ci_bt',ci_bt);



% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% estimating parameters
option = get(handles.popupmenu1,'Value');

switch option
    
case 1 

pars=evalin('base','pars');
minval=evalin('base','minval');
data=evalin('base','data');

% calculate confidence intervals
[ci_bt,se_bt]=model_ci_bt(pars,data);

assignin('base','ci_bt',ci_bt);
assignin('base','se_bt',se_bt);

% obtaining normal SE and CI
se=evalin('base','se');
ci=evalin('base','ci');

% Tabulate data
columns={'Parameters','SE','CI-low 95%','CI-high 95%','SE(BT)','CI(BT)-low 95%','CI(BT)-high 95%'};
table_data=horzcat(pars',se',ci(2,:)',ci(1,:)',se_bt',ci_bt(1,:)',ci_bt(2,:)');
position=[80/857 240/689 360/857 300/689];
set(handles.t,'Data',table_data,'ColumnName',columns,'Parent',handles.figure1,'Units','normalized','Position',position,'ColumnWidth','auto','FontUnits','normalized','FontSize',0.043);

% Enabling relevant buttons
set(handles.pushbutton4,'Enable','On');
set(handles.pushbutton5,'Enable','On');
set(handles.pushbutton6,'Enable','On');
set(handles.pushbutton7,'Enable','On');

case 2
            
     h=msgbox('SEs are computed only for the ML option','SE estimation','modal');
            
case 3
            
     h=msgbox('SEs are computed only for the ML option','SE estimation','modal');
            
end


% --------------------------------------------------------------------
function uipanel2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
