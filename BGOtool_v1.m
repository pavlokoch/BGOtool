function varargout = BGOtool_v1(varargin)
% BGOTOOL_V1 MATLAB code for BGOtool_v1.fig
%      BGOTOOL_V1, by itself, creates a new BGOTOOL_V1 or raises the existing
%      singleton*.
%
%      H = BGOTOOL_V1 returns the handle to a new BGOTOOL_V1 or the handle to
%      the existing singleton*.
%
%      BGOTOOL_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BGOTOOL_V1.M with the given input arguments.
%
%      BGOTOOL_V1('Property','Value',...) creates a new BGOTOOL_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BGOtool_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BGOtool_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BGOtool_v1

% Last Modified by GUIDE v2.5 08-Apr-2016 14:45:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BGOtool_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @BGOtool_v1_OutputFcn, ...
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


% --- Executes just before BGOtool_v1 is made visible.
function BGOtool_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BGOtool_v1 (see VARARGIN)

% get file list
handles.files = dir('*.dat');
set(handles.file_listbox,'String',{handles.files.name});
set(handles.info_textbox,'FontSize',8);

% exponential decay constant for BGO PMT
handles.tau = 0.3; %mus

% choose verbosity level
% 0 - report everything 
% 1 - exclude tables from reports

handles.verbosity = 1;
handles.min_energy_channel = 0;
handles.max_energy_channel = 4095;

% TGF windows
handles.HED_STW_D_1 = 300; %mus, short window 1
handles.HED_STW_D_2 = 1000; %mus, short window 2
handles.HED_STW_D_3 = 3000; %mus, short window 3
handles.HED_LTW_D = 25000; %mus, long window 

handles.BGO_Short_Threshold_1 = 1; % enable window 
handles.BGO_Short_Threshold_1 = 1; % enable window
handles.BGO_Short_Threshold_1 = 1; % enable window
handles.BGO_Short_Threshold_1 = 1; % enable window

% Vars for accepted counts algorithm
% Accepted Counts may be required to satisfy any of the following criteria: 
% - Have pulse height E between a lower limit E_L and an upper limit E_H, i.e. E_L ? E ? E_H
% - Have a time stamp T which, when compared with the time stamp of the previous count T_PREV, satisfies the relation T >= T_PREV + ?T
% Typical values for ?T are likely to be 0, 1 or 2 time stamp LSBs (nominally 0, 1 or 2 ?s).
handles.HED_STW_THR_L = 0; % 
handles.HED_STW_THR_U = 4095; % 
handles.HED_STW_AC_T = 0; % 255 disables the algorithm for short windows
handles.HED_LTW_THR_L = 0; % 
handles.HED_LTW_THR_U = 4095; % 
handles.HED_LTW_AC_T = 0; % 255 disables the algorithm for long window

% Choose default command line output for BGOtool_v1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%init context menu for lines, axes, etc.
init_contextmenu(handles);

% UIWAIT makes BGOtool_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BGOtool_v1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% basicaly main function

file_list = get(handles.file_listbox,'String');
selected_item = get(handles.file_listbox,'Value');

% clear all axes
arrayfun(@cla,findall(0,'type','axes'));
% clear info panel
set(handles.info_textbox,'String','');

% empty cell structure for collective data analysis
collective_data = {};

% loop over all selected files
for i = 1:numel(selected_item)

% create structures for messages, warnings and errors
handles.DAU.messages = {'DAU messages:'};
handles.DAU.warnings = {'DAU warnings:'};
handles.DAU.errors = {'DAU errors:'};

handles.DM0.messages = {'DM0 messages:'};
handles.DM0.warnings = {'DM0 warnings:'};
handles.DM0.errors = {'DM0 errors:'};

handles.DM1.messages = {'DM1 messages:'};
handles.DM1.warnings = {'DM1 warnings:'};
handles.DM1.errors = {'DM1 errors:'};

handles.DM2.messages = {'DM2 messages:'};
handles.DM2.warnings = {'DM2 warnings:'};
handles.DM2.errors = {'DM2 errors:'};  

   
% open, read and close file
    f = fopen(handles.files(selected_item(i)).name);
        if f,
            data = textscan(f,'%f %f %f %f %f');
            energy = data{:,1};
            flag = data{:,2};
            DMnum = data{:,3};  % DM number [0,1,2]
            time = data{:,4};
            ftime = data{:,5};
            fclose(f);
        else
            update_info(handles, 'error' ,'file not found');
        end
        
% flag meaning. ftime - fine time (fast time)
% 00 - 0 normal SCDP
% 01 - 1 fast
% 10 - 2 overflow 
% 10 and ftime MSB - valley
% 11 - 3 ADC sample             

% --------- process and plot -------------

% Preliminary analysis
    handles.DAU.messages{end+1} = ['Analysis of ',handles.files(selected_item(i)).name];

% energy range 0 - 4096
if sum(energy < 0 | energy > 4095)
     handles.DAU.warnings{end+1} = 'Energy out of range (0-4095).';
end

% check for zero energy
if sum(energy==0)
    handles.DAU.warnings{end+1} = ['Zero energy detected in ', num2str(sum(energy==0),'%3.0f'), ' SCDPs.'];
end

% fraction of Fast/OVF events
    handles.DAU.messages{end+1} = ['Fraction of Fast/Valley events ', num2str(100*sum(flag==1)/numel(flag),'%2.2f'),'%'];

% check for continuous operation mode
if numel(find(flag == 3)),
    handles.DAU.messages{end+1} = 'ADC sample detected (flag==3)! ';    
end

% find unique channels
uniq_DMnum = unique(DMnum);

% report
handles.DAU.messages{end+1} = ['DMs found: ',num2str(uniq_DMnum')];

% check if number of channels equals to three
if numel(uniq_DMnum) ~= 3, 
    handles.DAU.warnings{end+1} = 'Missing DMs';
    set(handles.dm0_checkbox,'Value',sum(uniq_DMnum==0));
    set(handles.dm1_checkbox,'Value',sum(uniq_DMnum==1));
    set(handles.dm2_checkbox,'Value',sum(uniq_DMnum==2));
    % Update handles structure
    guidata(hObject, handles);
end; 

%
% All other checks must go here, or later in another function
%

% convert cell to numbers
data = cell2mat(data);

% process fast valley events in each DM
% function substracts all valleys from their fasts and remove valleys from
% dataset
[data_DM0, handles.DM0] = fast_valley_processing(handles,handles.DM0,data(DMnum == 0,:));
[data_DM1, handles.DM1] = fast_valley_processing(handles,handles.DM1,data(DMnum == 1,:));
[data_DM2, handles.DM2] = fast_valley_processing(handles,handles.DM2,data(DMnum == 2,:));

% *************************************************************************
% remove fast SCDP from analysis if enabled
if ~get(handles.fast_valley_checkbox,'Value')
data_DM0 = data_DM0(data_DM0(:,2) == 0,:);
data_DM1 = data_DM1(data_DM1(:,2) == 0,:);
data_DM2 = data_DM2(data_DM2(:,2) == 0,:);
handles.DAU.warnings{end+1} = 'ALL FAST EVENTS EXCLUDED!';
end
% *************************************************************************

% to show/not to show (checkbox on Control panel)
show_DM0 = get(handles.dm0_checkbox,'Value');
show_DM1 = get(handles.dm1_checkbox,'Value');
show_DM2 = get(handles.dm2_checkbox,'Value');

% if anything to show 
if show_DM0 || show_DM1 || show_DM2

    if show_DM0,
        % make linear time scale
        data_DM0(:,4) = linear_time_scale(data_DM0(:,4),data_DM0(:,5));
        % remove ftime from further analysis
        data_DM0(:,5) = [];
        % if individual analysis required (checkbox on Options panel)
        if ~get(handles.collective_checkbox,'Value')
            [status, handles.DM0] = process(handles, handles.DM0, data_DM0, 'r');
        else
        % else collect data for further analysis after this loop
        collective_data(end+1) = {data_DM0}; 
        end
    end 
    if show_DM1,
        % make linear time scale
        data_DM1(:,4) = linear_time_scale(data_DM1(:,4),data_DM1(:,5));
        % remove ftime from further analysis
        data_DM1(:,5) = [];
        % if individual analysis required (checkbox on Options panel)
        if ~get(handles.collective_checkbox,'Value')
            [status, handles.DM1] = process(handles, handles.DM1, data_DM1, 'g');
        else
        % else collect data for further analysis after this loop
        collective_data(end+1) = {data_DM1}; 
        end 
    end
    if show_DM2,
        % make linear time scale
        data_DM2(:,4) = linear_time_scale(data_DM2(:,4),data_DM2(:,5));
        % remove ftime from further analysis
        data_DM2(:,5) = [];
        % if individual analysis required (checkbox on Options panel)
        if ~get(handles.collective_checkbox,'Value')
            [status, handles.DM2] = process(handles, handles.DM2, data_DM2, 'b');
        else
        % else collect data for further analysis after this loop
        collective_data(end+1) = {data_DM2}; 
        end
    end

%if nothing to show clear axes
else
   arrayfun(@cla,findall(0,'type','axes')); 
end    

% print all messages to info panel
report(handles);

end     % end loop over all selected files

% if collective analysis enabled on Options panel 
if get(handles.collective_checkbox,'Value')
           % stack together 
           collective_data = vertcat(collective_data{1:end});
           % then sort rows in time 
           collective_data = sortrows(collective_data,4);
           % and fimally procees
           handles.DAU.messages{end+2} = '***** Collective data analysis: *****';
           [status, handles.DAU] = process(handles, handles.DAU, collective_data, 'b');
           report(handles);
end

function [status, DM] = process(handles, DM, data, color)
% main function to process and display all data
% DM - communication pipeline 

% split
energy = data(:,1);
flag = data(:,2);
time = data(:,4);

% get intervals in mus
intervals_DM = diff(time);

% report total observation time
    DM.messages{end+1} = ['Total observation time: ',num2str(range(time)/1e6,'%3.0f'),' sec'];

% check if time monotonously increases
if sum(intervals_DM < 0),
    DM.errors{end+1} = 'TIME IS NOT MONOTONOUSLY INCREASING!';
    DM.errors{end+1} = 'Check following times:';
    DM.errors{end+1} = num2str(time(intervals_DM<0),'%9.0f');
    line(handles.data_axes,[time(intervals_DM<0),time(intervals_DM<0)],[0,handles.max_energy_channel])
end

% ---oooOOO Timing OOOooo---
% prepare histogram
time_bin = 1000; % in mus
% define edges
edges = 0:time_bin:max(intervals_DM);
% check if axes is hold
if ~ishold(handles.timing_axes),hold(handles.timing_axes);end;
% make histogram
counts = histcounts(intervals_DM,edges); 
% normalize
counts = counts/sum(counts);

    [l,gof] = plot_time_intervals(handles, edges(1:end-1)+time_bin/2, counts, time_bin, color);
    DM.messages{end+1} = ['Arrival rate (lambda): ',num2str(l*1e3,'%3.2f'),' counts/ms'];
    DM.messages{end+1} = ['Poisson gof (rsquare): ',num2str(gof)];

% ---oooOOO Deposited spectra OOOooo---
% prepare histogram
energy_bin = 100; % in channels    
% define edges    
edges = handles.min_energy_channel:energy_bin:handles.max_energy_channel;
% make histogram
counts = histcounts(energy,edges);
% normalize
counts = counts/sum(counts)/energy_bin;
% check if axes is hold
if ~ishold(handles.spectra_axes),hold(handles.spectra_axes);end;

    plot_energy_spectra(handles, edges(1:end-1)+energy_bin/2, counts, energy_bin, range(time)/1e6, color)

% ---oooOOO Primary spectra OOOooo---
% only for test purposes. incident angel needs to be provided for correct
% DRM usage. 
% All DRMs binned from 0.01 to 100 MeV in 40 log bins. 
energy_edges = [0.01 0.0125893 0.0158489 0.0199526 0.0251189 0.0316228 0.0398107 ... 
         0.0501187 0.0630957 0.0794328 0.1 0.125893 0.158489 0.199526 0.251189 ...
         0.316228 0.398107 0.501187 0.630957 0.794328 1 1.25893 1.58489 1.99526 ...
         2.51189 3.16228 3.98107 5.01187 6.30957 7.94328 10 12.5893 15.8489 ... 
         19.9526 25.1189 31.6228 39.8107 50.1187 63.0957 79.4328 100]; % in MeV

% scale to channels (ref ASIM-UB-UBINT-RP-006 BGO DM spectroscopic tests.pdf page 13)
% for a moment lets assume 0.01 - min energy channel and 100 hits max
% energy channel
ch_edges = energy_edges*(handles.max_energy_channel - handles.min_energy_channel)/max(energy_edges); 
% create histogram
counts = histcounts(energy,ch_edges);
% load test DRM for 0 0 angles. Correct DRM should be loaded according to the incident
% angles.
load('test_drm.mat');
% check if axes is hold
if ~ishold(handles.primary_axes),hold(handles.primary_axes);end;
% convert deposited to primary using DRM
pr_spectra = deposited2primary(drm,counts);

    stairs(handles.primary_axes, energy_edges(1:end-1), pr_spectra, color);

set(handles.primary_axes,'XScale','log');    
xlabel(handles.primary_axes, 'Energy [MeV]');
ylabel(handles.primary_axes, 'Counts / channel / sec / cm^2');
legend(handles.primary_axes,'Primary spectrum');

% ---oooOOO Real-time data OOOooo---
% check if axes is hold
if ~ishold(handles.data_axes),hold(handles.data_axes);end;

    ln = stem(handles.data_axes, time, energy, ['.',color]);
    % plot fast events if any
    stem(handles.data_axes, time(flag == 1), energy(flag == 1), ['*',color]);
    set(ln, 'uicontextmenu',handles.hcmenu_timing); 

xlabel(handles.data_axes, 'Time [\mus]');
ylabel(handles.data_axes, 'Energy [ch]');


% no errors
status = 1;


% % --- Executes on button press in dm0_checkbox.
function dm0_checkbox_Callback(hObject, eventdata, handles)
% pause(0.7);
% update_button_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in dm1_checkbox.
function dm1_checkbox_Callback(hObject, eventdata, handles)
% pause(0.7);
% update_button_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in dm2_checkbox.
function dm2_checkbox_Callback(hObject, eventdata, handles)
% pause(0.7);
% update_button_Callback(hObject, eventdata, handles)

% --- Executes on selection change in file_listbox.
function file_listbox_Callback(hObject, eventdata, handles)
% clear all axes
arrayfun(@cla,findall(0,'type','axes'));
% update
update_button_Callback(hObject, eventdata, handles)


function init_contextmenu(handles)
% general context menu
handles.hcmenu = uicontextmenu;
% for timing plots
handles.hcmenu_timing = uicontextmenu;
% for info pane;
handles.hcmenu_info = uicontextmenu;


% Update handles structure
guidata(gcf, handles);
%export data to MATLAB workspace
uimenu(handles.hcmenu, 'Label', 'Export to WS', 'Callback', {@export_data,handles});
uimenu(handles.hcmenu_timing, 'Label', 'Export to WS', 'Callback', {@export_data,handles});

uimenu(handles.hcmenu_timing, 'Label', 'Fit', 'Callback', {@exp_decay,handles});

hc_info = uimenu(handles.hcmenu_info, 'Label', 'verbosity');
uimenu('Parent',hc_info, 'Label','0','Callback',{@set_verbosity,handles});
uimenu('Parent',hc_info, 'Label','1','Callback',{@set_verbosity,handles});
uimenu('Parent',hc_info, 'Label','2','Callback',{@set_verbosity,handles});
set(handles.info_textbox, 'uicontextmenu', handles.hcmenu_info);

function set_verbosity(hObject,eventdata,handles)
handles.verbosity = str2double(eventdata.Source.Label);
% Update handles structure
guidata(gcf, handles);
update_button_Callback(hObject, eventdata, handles)

function export_data(~,~,handles)
assignin('base', 'data', [get(gco,'XData');get(gco,'YData')]');

function exp_decay(~, ~, handles)
point = ginput(1);
x = 0:0.1:5;
% time decay constant 0.3 mus
plot(gca,point(1)+x,point(2)*exp(-x/handles.tau),'r');

function [lambda,gof] = plot_time_intervals(handles, edges, counts, time_bin, style)
    % plot measured interval distribution
    l = stairs(handles.timing_axes,edges,counts,style);
    % assign context menu
    set(l, 'uicontextmenu',handles.hcmenu); 
    % exponential fit
    [ft,gof] = fit(edges',counts',[num2str(time_bin),'*lambda*exp(-lambda*x)'],'StartPoint',0);
    % focus on axes
    axes(handles.timing_axes);
    %plot fit
    plot(ft,['--',style]); legend('off');
    lambda = ft.lambda;
    gof = gof.rsquare;
    
% scale and labels

set(handles.timing_axes,'XScale','log');
set(handles.timing_axes,'YScale','log');

xlabel(handles.timing_axes, 'Time intervals [\mus]');
ylabel(handles.timing_axes, 'Counts');

    
function plot_energy_spectra(handles, edges, counts, energy_bin, observation_time, style)
    % plot measured spectra
    l = stairs(handles.spectra_axes,edges,counts./observation_time,style);
    % assign context menu    
    set(l, 'uicontextmenu',handles.hcmenu);
    
set(handles.spectra_axes,'YScale','log');

xlabel(handles.spectra_axes, 'Channel');
ylabel(handles.spectra_axes, 'Counts / channel / sec');

legend(handles.spectra_axes, 'Deposited spectrum');
    
function [processed_data, DM] = fast_valley_processing(handles,DM, data)
% Find complete pairs of fast/valley events to correct time, energy and
% swap them. Pairs must comply:
% - first row is the fast event (flag==1)
% - second row is the valley event (flag==2 and ftime MSB==1)
% - absolute time of the second event must be lower than of the first one

% first find all fasts (flag==1)
ind = find(data(:,2) == 1);
fasts = data(ind,:);
% valleys should follow fasts, be within 2 mus and have MSB == 1. Otherwise
% reject.
valleys = zeros(size(fasts));
incorrect_valleys = false([numel(ind),1]);
for i = 1:numel(ind)
    temp = data(ind(i)+1,:);
    valleys(i,:) = temp;
    if (temp(1,4) - fasts(i,4) <= 2) && temp(1,5) >= 128,
        incorrect_valleys(i) = 0;
    else
        incorrect_valleys(i) = 1;
    end
end

% check if all fast events have valley pair
if sum(incorrect_valleys) && handles.verbosity == 0,
    DM.warnings{end+1} = 'incorrect valleys detected:';
    DM.warnings{end+1} = num2str(valleys(incorrect_valleys,:));
end

% correct valley time
valleys(:,5) = valleys(:,5) - 128;
% check valleys time
incorrect_time = fasts(:,4) + fasts(:,5)*027.77777e-3 < valleys(:,4) + valleys(:,5)*027.77777e-3;
if sum(incorrect_time) && handles.verbosity == 0,
   DM.warnings{end+1} = 'incorrect valley time:'; 
   DM.warnings{end+1} = num2str(valleys(incorrect_time,:));
end

% substract valley from its fast if enabled
% Pulse decays exponentially with time constant approx 300 ns at 20 deg [DOI:10.1109/I2MTC.2012.6229186]
if get(handles.tail_comp_checkbox,'Value')
fasts(:,1) = round(fasts(:,1) - valleys(:,1).*exp(-(fasts(:,4) - valleys(:,4) + 027.77777e-3*(fasts(:,5) - valleys(:,5)))/handles.tau));
end
% check if energy <= 0
incorrect_energy = fasts(:,1) <= 0;
if sum(incorrect_energy) && handles.verbosity == 0,
   DM.warnings{end+1} = 'incorrect Fast energy:';
   DM.warnings{end+1} = num2str(round(fasts(incorrect_energy,:)));
end
% remove all incorrect fasts
fasts(incorrect_time,:) = NaN;
fasts(incorrect_energy,:) = NaN;
fasts(incorrect_valleys,:) = NaN;

processed_data = data;
% replace with corrected fast events
processed_data(ind,:) = fasts;
% remove valleys from the data
processed_data(ind+1,:) = [];

if sum(isnan(processed_data(:,1))),
    DM.warnings{end+1} = [num2str(sum(isnan(fasts(:,1)))),' Fast SCDPs removed'];
end


function report(handles)  
% show all messages from DAUs and DMs
current_text = {get(handles.info_textbox,'String')};
if numel(current_text{1}) > 1,
    current_text = current_text{:};
end
current_text(end+1:end+numel(handles.DAU.messages)) = handles.DAU.messages;
current_text(end+1:end+numel(handles.DAU.warnings)) = handles.DAU.warnings;
current_text(end+1:end+numel(handles.DAU.errors)) = handles.DAU.errors;
current_text(end+1) = {'- - - - - - - - - - - - -'};

if get(handles.dm0_checkbox,'Value')
current_text(end+1:end+numel(handles.DM0.messages)) = handles.DM0.messages;
current_text(end+1:end+numel(handles.DM0.warnings)) = handles.DM0.warnings;
current_text(end+1:end+numel(handles.DM0.errors)) = handles.DM0.errors;
current_text(end+1) = {'- - - - - - - - - - - - -'};
end

if get(handles.dm1_checkbox,'Value')
current_text(end+1:end+numel(handles.DM1.messages)) = handles.DM1.messages;
current_text(end+1:end+numel(handles.DM1.warnings)) = handles.DM1.warnings;
current_text(end+1:end+numel(handles.DM1.errors)) = handles.DM1.errors;
current_text(end+1) = {'- - - - - - - - - - - - -'};
end

if get(handles.dm2_checkbox,'Value')
current_text(end+1:end+numel(handles.DM2.messages)) = handles.DM2.messages;
current_text(end+1:end+numel(handles.DM2.warnings)) = handles.DM2.warnings;
current_text(end+1:end+numel(handles.DM2.errors)) = handles.DM2.errors;
current_text(end+1) = {'- - - - - - - - - - - - -'};
end

set(handles.info_textbox,'String',current_text);




