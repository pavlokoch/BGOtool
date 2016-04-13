function varargout = triggers(varargin)
% triggers MATLAB code for triggers.fig
%
%      H = triggers returns the handle to a new triggers or the handle to
%      the existing singleton*.
%
%      triggers('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in triggers.M with the given input arguments.
%
%      triggers('Property','Value',...) creates a new triggers or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before triggers_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to triggers_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help triggers

% Last Modified by GUIDE v2.5 12-Apr-2016 11:53:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @triggers_OpeningFcn, ...
                   'gui_OutputFcn',  @triggers_OutputFcn, ...
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


% --- Executes just before triggers is made visible.
function triggers_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to triggers (see VARARGIN)

% get file list
handles.files = dir('*.dat');
set(handles.file_listbox,'String',{handles.files.name});
set(handles.info_textbox,'FontSize',8);

% exponential decay constant for BGO PMT
handles.tau = 0.3; %mus

% choose verbosity level
% 0 - report everything 
% 1 - exclude tables from reports
% 2 - only triggers
handles.verbosity = 2;

handles.min_energy_channel = 0;
handles.max_energy_channel = 4095;

% TGF windows
handles.HED_STW_D_1 = 300; %mus, short window 1
handles.HED_STW_D_2 = 1000; %mus, short window 2
handles.HED_STW_D_3 = 3000; %mus, short window 3
handles.HED_LTW_D = 25000; %mus, long window 

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

% Choose default command line output for triggers
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%init context menu for lines, axes, etc.
init_contextmenu(handles);

% UIWAIT makes triggers wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = triggers_OutputFcn(hObject, eventdata, handles) 
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
if handles.verbosity < 2
    handles.DAU.messages{end+1} = ['Analysis of ',handles.files(selected_item(i)).name];
end

% energy range 0 - 4096
if sum(energy < 0 | energy > 4095) && handles.verbosity < 2
     handles.DAU.warnings{end+1} = 'Energy out of range (0-4095).';
end

% check for zero energy
if sum(energy==0) && handles.verbosity < 2
    handles.DAU.warnings{end+1} = ['Zero energy detected in ', num2str(sum(energy==0),'%3.0f'), ' SCDPs.'];
end

% fraction of Fast/OVF events
if  handles.verbosity < 2
    handles.DAU.messages{end+1} = ['Fraction of Fast/Valley events ', num2str(100*sum(flag==1)/numel(flag),'%2.2f'),'%'];
end
    
% check for continuous operation mode
if numel(find(flag == 3)) && handles.verbosity < 2,
    handles.DAU.messages{end+1} = 'ADC sample detected (flag==3)! ';    
end

% find unique channels
uniq_DMnum = unique(DMnum);

% report
if  handles.verbosity < 2,
handles.DAU.messages{end+1} = ['DMs found: ',num2str(uniq_DMnum')];
end

% check if number of channels equals to three
if numel(uniq_DMnum) ~= 3 && handles.verbosity < 2, 
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

end     % end loop over all selected files

% if collective analysis enabled on Options panel 
if get(handles.collective_checkbox,'Value')
           % stack together 
           collective_data = vertcat(collective_data{1:end});
           % then sort rows in time 
           collective_data = sortrows(collective_data,4);
           % and fimally procees
           handles.DAU.messages{end+2} = '***** Collective data analysis: *****';
           handles.DAU.warnings{end+1} = 'Please check each DM individually for time errors!';
           [status, handles.DAU] = process(handles, handles.DAU, collective_data, 'b');
           report(handles);
else
% just report
report(handles);
end

function [status, DM] = process(handles, DM, data, color)
% main function to process and display all data
% DM - communication pipeline 

% split
energy = data(:,1);
flag = data(:,2);
time = data(:,4);

dt = diff(time);
if sum(dt < 0),
    DM.errors{end+1} = 'TIME IS NOT MONOTONOUSLY INCREASING!';
    DM.errors{end+1} = 'Check following times:';
    DM.errors{end+1} = num2str(time(dt<0),'%9.0f');
    line(handles.data_axes,[time(dt<0),time(dt<0)],[0,handles.max_energy_channel])
end

% total observation time
    DM.messages{end+1} = ['Total observation time: ',num2str(round(range(time))/1e6),' sec'];


% ---oooOOO Real-time data OOOooo---
% check if axes is hold
if ~ishold(handles.data_axes),hold(handles.data_axes);end;

    ln = stem(handles.data_axes, time, energy, ['.',color]);
    % plot fast events if any
    stem(handles.data_axes, time(flag == 1), energy(flag == 1), ['*',color]);
    set(ln, 'uicontextmenu',handles.hcmenu_timing); 

xlabel(handles.data_axes, 'Time [\mus]');
ylabel(handles.data_axes, 'Energy [ch]');

% ---oooOOO TGF triggers OOOooo---
% Trigger algorithm strictly accordingly to documentation.
% Trigger flags down
handles.BGO_Trigger_Flag_Short_1 = false([1,numel(time)]);
handles.BGO_Trigger_Flag_Short_2 = false([1,numel(time)]);
handles.BGO_Trigger_Flag_Short_3 = false([1,numel(time)]);
handles.BGO_Trigger_Flag_Long = false([1,numel(time)]);

% variables for in-flight background rate estimation
% two backgrounds now implemented - for short and for long timewindows
    TCP_counter_short = 0;
    TCP_position_short = 0;
    TCP_counter_long = 0;
    TCP_position_long = 0;
    buffSize = 8;
    C_CircBuff_short = nan(1,buffSize);
    C_CircBuff_long = nan(1,buffSize);

    % number of selected modules to compensate for bg estimation
    % if not complete mxgs data selected then extrapolate background rate
    dm_num = numel(get(handles.file_listbox,'Value')) * sum(get(handles.dm0_checkbox,'Value') + ...
                                                            get(handles.dm1_checkbox,'Value') + ...
                                                            get(handles.dm2_checkbox,'Value'));
    % for beginning let background be average count rate
    bg_short = numel(time)/range(time)*1e3/dm_num*12; %counts/ms    
    bg_long = numel(time)/range(time)*1e3/dm_num*12; %counts/ms    

                                                        
% include first record
BGO_S_TS_Buffer(1) = time(1);
BGO_L_TS_Buffer(1) = time(1);
% loop over all entries starting from second element
for i = 2:numel(energy)
   
    % only start if background is estimated
    if bg_short > 0 || bg_long > 0
   
    % first short timewindow processing 
    % check energy range
   if energy(i) >= handles.HED_STW_THR_L && energy(i) <= handles.HED_STW_THR_U 
       % if short window AC algorithm enabled
       if handles.HED_STW_AC_T < 255 
           if (time(i) - time(i-1) >= handles.HED_STW_AC_T)
               BGO_S_TS_Buffer(end+1) = time(i);     
           end
       else
           BGO_S_TS_Buffer(end+1) = time(i);
       end
       % check short window 1 for trigger (if enabled and enough elements in buffer)
       if ~handles.HED_STW_D_1 == 0 && (numel(BGO_S_TS_Buffer) > trigger_lookups(bg_short, 1))
          T_thresh = BGO_S_TS_Buffer(end - trigger_lookups(bg_short, 1));
          T_delta = BGO_S_TS_Buffer(end) - T_thresh;
            if T_delta <= handles.HED_STW_D_1, handles.BGO_Trigger_Flag_Short_1(i) = true;
            end
       end
       % check short window 2 for trigger (if enabled and enough elements in buffer)
       if ~handles.HED_STW_D_2 == 0 && (numel(BGO_S_TS_Buffer) > trigger_lookups(bg_short, 2))
          T_thresh = BGO_S_TS_Buffer(end - trigger_lookups(bg_short, 1));
          T_delta = BGO_S_TS_Buffer(end) - T_thresh;
            if T_delta <= handles.HED_STW_D_2, handles.BGO_Trigger_Flag_Short_2(i) = true;
            end
       end
       % check short window 3 for trigger (if enabled and enough elements in buffer)
       if ~handles.HED_STW_D_3 == 0 && (numel(BGO_S_TS_Buffer) > trigger_lookups(bg_short, 3))
          T_thresh = BGO_S_TS_Buffer(end - trigger_lookups(bg_short, 1));
          T_delta = BGO_S_TS_Buffer(end) - T_thresh;
            if T_delta <= handles.HED_STW_D_3, handles.BGO_Trigger_Flag_Short_3(i) = true;
            end
       end
       
   end
   % long timewindow processing     
   if energy(i) >= handles.HED_LTW_THR_L && energy(i) <= handles.HED_LTW_THR_U
       % if long window AC algorithm enabled
       if handles.HED_LTW_AC_T < 255 
           if (time(i) - time(i-1) >= handles.HED_LTW_AC_T)
               BGO_L_TS_Buffer(end+1) = time(i);
           end
       else
           BGO_L_TS_Buffer(end+1) = time(i);
       end
      % check long window for trigger (if enabled and enough elements in buffer)
       if ~handles.HED_LTW_D == 0 && (numel(BGO_L_TS_Buffer) > trigger_lookups(bg_long, 4))
          T_thresh = BGO_L_TS_Buffer(end - trigger_lookups(bg_long, 1));
          T_delta = BGO_L_TS_Buffer(end) - T_thresh;
            if T_delta <= handles.HED_LTW_D, handles.BGO_Trigger_Flag_Long(i) = true;
            end
       end
   end 
   
   % background ratemeter for short window
   % if TCP arrived
   if floor(BGO_S_TS_Buffer(end) / (1000000*(TCP_counter_short + 1))) > 0
       TCP_counter_short = TCP_counter_short + 1;
       C_CircBuff_short = [C_CircBuff_short(2:end) (numel(BGO_S_TS_Buffer) - TCP_position_short)]; % background per second
       TCP_position_short = numel(BGO_S_TS_Buffer);  % save last TCP position
       if TCP_counter_short >= 8  % if 8 seconds passed
          % update background for the next second
          % see formula in Data Acq Algorithm
          S = sum(C_CircBuff_short);
          ST = -sum(C_CircBuff_short(1:end-1).*[7:-1:1]);
          bg_short = (42*S + 9*ST)/84/1e3;
          % compensate for number of selected modules
          bg_short = bg_short/dm_num*12;
       end
   end
   % background ratemeter for long window
   % if TCP arrived
   if floor(BGO_L_TS_Buffer(end) / (1000000*(TCP_counter_long + 1))) > 0
       TCP_counter_long = TCP_counter_long + 1;
       C_CircBuff_long = [C_CircBuff_long(2:end) (numel(BGO_L_TS_Buffer) - TCP_position_long)]; % background per second
       TCP_position_long = numel(BGO_L_TS_Buffer);  % save last TCP position
       if TCP_counter_long >= 8  % if 8 seconds passed
          % update background for the next second
          % see formula in Data Acq Algorithm
          S = sum(C_CircBuff_long);
          ST = -sum(C_CircBuff_long(1:end-1).*[7:-1:1]);
          bg_long = (42*S + 9*ST)/84/1e3;
          % compensate for number of selected modules
          bg_long = bg_long/dm_num*12;
       end
   end
   
   
   % CZT-BGO Coincidence comes here
   % ...
    end % end if bg > 0
end

% check if axes is hold
if ~ishold(handles.SW1_axes),hold(handles.SW1_axes);end;
if ~ishold(handles.SW2_axes),hold(handles.SW2_axes);end;
if ~ishold(handles.SW3_axes),hold(handles.SW3_axes);end;
if ~ishold(handles.LW_axes),hold(handles.LW_axes);end;

% plot triggerst
stem(handles.SW1_axes, time, handles.BGO_Trigger_Flag_Short_1,['.',color])
stem(handles.SW2_axes, time, handles.BGO_Trigger_Flag_Short_2,['.',color])
stem(handles.SW3_axes, time, handles.BGO_Trigger_Flag_Short_3,['.',color])
stem(handles.LW_axes, time, handles.BGO_Trigger_Flag_Long,['.',color])

% axes limits
xlim(handles.data_axes, [0, max(time)]);
ylim(handles.SW1_axes, [0,2]);
ylim(handles.SW2_axes, [0,2]);
ylim(handles.SW3_axes, [0,2]);
ylim(handles.LW_axes, [0,2]);

% axes labels
ylabel(handles.SW1_axes, 'SW1');
ylabel(handles.SW2_axes, 'SW2');
ylabel(handles.SW3_axes, 'SW3');
ylabel(handles.LW_axes, 'LW');

% link axes together
linkaxes([handles.data_axes,handles.SW1_axes,handles.SW2_axes,handles.SW3_axes,handles.LW_axes],'x');

% report
DM.messages{end+1} = ['Short Window 1 triggers (total): ', num2str(sum(handles.BGO_Trigger_Flag_Short_1))];
DM.messages{end+1} = ['Short Window 1 triggers (per hour): ', num2str(round(sum(handles.BGO_Trigger_Flag_Short_1)/range(time)*1e6*3600))];
DM.messages{end+1} = ['Short Window 2 triggers (total): ', num2str(sum(handles.BGO_Trigger_Flag_Short_2))];
DM.messages{end+1} = ['Short Window 2 triggers (per hour): ', num2str(round(sum(handles.BGO_Trigger_Flag_Short_2)/range(time)*1e6*3600))];
DM.messages{end+1} = ['Short Window 3 triggers (total): ', num2str(sum(handles.BGO_Trigger_Flag_Short_3))];
DM.messages{end+1} = ['Short Window 3 triggers (per hour): ', num2str(round(sum(handles.BGO_Trigger_Flag_Short_3)/range(time)*1e6*3600))];
DM.messages{end+1} = ['Long Window triggers (total): ', num2str(sum(handles.BGO_Trigger_Flag_Long))];
DM.messages{end+1} = ['Long Window triggers (per hour): ', num2str(round(sum(handles.BGO_Trigger_Flag_Long)/range(time)*1e6*3600))];

% no errors
status = 1;

% --- Executes on selection change in file_listbox.
function file_listbox_Callback(hObject, eventdata, handles)
% clear all axes
arrayfun(@cla,findall(0,'type','axes'));
% clear info panel
set(handles.info_textbox,'String','');


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

function export_data(~,~,handles)
assignin('base', 'data', [get(gco,'XData');get(gco,'YData')]');

function exp_decay(~, ~, handles)
point = ginput(1);
x = 0:0.1:5;
% time decay constant 0.3 mus
plot(gca,point(1)+x,point(2)*exp(-x/handles.tau),'r');

    
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
% fasts(incorrect_time,:) = NaN;
% fasts(incorrect_energy,:) = NaN;
% fasts(incorrect_valleys,:) = NaN;

processed_data = data;
% replace with corrected fast events
processed_data(ind,:) = fasts;
% remove valleys from the data
processed_data(ind+1,:) = [];

if sum(isnan(processed_data(:,1)))  && handles.verbosity < 2,
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



function SW1_edit_Callback(hObject, eventdata, handles)
handles.HED_STW_D_1 = str2double(get(handles.SW1_edit,'String')); %mus, short window 1 duration
% Update handles structure
guidata(gcf, handles);

function SW2_edit_Callback(hObject, eventdata, handles)
handles.HED_STW_D_2 = str2double(get(handles.SW2_edit,'String')); %mus, short window 2 duration
% Update handles structure
guidata(gcf, handles);

function SW3_edit_Callback(hObject, eventdata, handles)
handles.HED_STW_D_3 = str2double(get(handles.SW3_edit,'String')); %mus, short window 3 duration
% Update handles structure
guidata(gcf, handles);

function LW_edit_Callback(hObject, eventdata, handles)
handles.HED_LTW_D = str2double(get(handles.LW_edit,'String')); %mus, long window duration
% Update handles structure
guidata(gcf, handles);


function STW_THR_L_edit_Callback(hObject, eventdata, handles)
handles.HED_STW_THR_L = str2double(get(handles.STW_THR_L_edit,'String')); % short time window low energy threshold
% Update handles structure
guidata(gcf, handles);

function STW_THR_U_edit_Callback(hObject, eventdata, handles)
handles.HED_STW_THR_U = str2double(get(handles.STW_THR_U_edit,'String')); % short time window up energy threshold
% Update handles structure
guidata(gcf, handles);

function STW_AC_T_edit_Callback(hObject, eventdata, handles)
handles.HED_STW_AC_T = str2double(get(handles.STW_AC_T_edit,'String')); % short time window accepeted count time
% Update handles structure
guidata(gcf, handles);

function LTW_THR_L_edit_Callback(hObject, eventdata, handles)
handles.HED_LTW_THR_L = str2double(get(handles.LTW_THR_L_edit,'String')); % long time window low energy threshold
% Update handles structure
guidata(gcf, handles);

function LTW_THR_U_edit_Callback(hObject, eventdata, handles)
handles.HED_LTW_THR_U = str2double(get(handles.LTW_THR_U_edit,'String')); % long time window up energy threshold
% Update handles structure
guidata(gcf, handles);

function LTW_AC_T_edit_Callback(hObject, eventdata, handles)
handles.HED_LTW_AC_T = str2double(get(handles.LTW_AC_T_edit,'String')); % long time window accepeted count time
% Update handles structure
guidata(gcf, handles);
