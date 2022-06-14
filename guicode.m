function varargout = gui(varargin)
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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

function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

function pushbutton1_Callback(hObject, eventdata, handles)
[filename filepath]=uigetfile({'*.mat'},'choose your file');
set(handles.edit1,'string',[filepath filename]) %show in edit1
road=get(handles.edit1,'String');
handles.road=road;
guidata(hObject, handles)

function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton2_Callback(hObject, eventdata, handles)
road=handles.road;
a = load(road);  % load the file road

channel_num = str2double(get(handles.edit2,'string'));
handles.channel_num = channel_num;
guidata(hObject,handles)

% get the filename in a
b = fieldnames(a);
name = strcat('a.',b(1,1));
s = eval(name{1});
handles.s = s;
guidata(hObject,handles)

% plot the signals
axes(handles.axes1);
maxsample = max(max(s))*1.1;
handles.maxsample = maxsample;
guidata(hObject,handles)
for i=1:channel_num 
    plot(s(:,i) - maxsample*(i-1)*1,'LineWidth',.6)
    hold on
end
yticks([-(maxsample*(channel_num-1)):maxsample:0])
for i=1:channel_num
    ticklabels{i} = strcat('ch',num2str(channel_num-i+1));
end
yticklabels(ticklabels)
xlim([-inf,inf])
ylim([-(maxsample*channel_num),maxsample])
ylabel('Channels')
xlabel('Samples')

function pushbutton3_Callback(hObject, eventdata, handles)
% clean the axes2
try
    delete(allchild(handles.axes2));
end

channel_num = handles.channel_num;
s = handles.s;
maxsample = handles.maxsample;
fs = str2double(get(handles.edit7,'string'));
handles.fs = fs;
guidata(hObject, handles)

% plot the preprocessing signals
axes(handles.axes2);  
% Butterworth filter
Wc = 2 * 1 / fs;  % 1Hz
[b, a] = butter(2, Wc);  % second-order
for i=1:channel_num 
    y = s(:,i);
    y_Filter = filter(b, a, y);
    bufi(:,i) = y_Filter;
    plot(y_Filter - maxsample*(i-1)*1,'LineWidth',.6)
    hold on
end
yticks([-(maxsample*(channel_num-1)):maxsample:0])
for i=1:channel_num
    ticklabels{i} = strcat('ch',num2str(channel_num-i+1));
end
yticklabels(ticklabels)
xlim([-inf,inf])
ylim([-(maxsample*channel_num),maxsample])
ylabel('Channels')
xlabel('Samples')
handles.bufi = bufi;
guidata(hObject, handles)
now = bufi;
handles.now = now;
guidata(hObject, handles)

function pushbutton4_Callback(hObject, eventdata, handles)
try
    delete(allchild(handles.axes2));
end

% plot the preprocessing signals
axes(handles.axes2);
bufi = handles.bufi;
fs = handles.fs;
channel_num = handles.channel_num;
maxsample = handles.maxsample;

gesture_num = str2double(get(handles.edit3,'string'));
handles.gesture_num = gesture_num;
guidata(hObject, handles)
repeat = str2double(get(handles.edit4,'string'));
handles.repeat = repeat;
guidata(hObject, handles)
action_time = str2double(get(handles.edit5,'string'));
handles.action_time = action_time;
guidata(hObject, handles)
rest_time = str2double(get(handles.edit6,'string'));
handles.rest_time = rest_time;
guidata(hObject, handles)

aftercut = action_time*fs*0.6; % remain the middle 60 percent of the signals
handles.aftercut = aftercut;
guidata(hObject, handles)
beforecut = (action_time+rest_time)*fs;
start = 1+action_time*fs/2-aftercut/2;
over = action_time*fs/2+aftercut/2;

for i=1:channel_num
    for k=1:gesture_num
        for j = 1:repeat
            % n = j+(k-1)*repeat,this is nth action in total
            n = j+(k-1)*repeat;
            cut(1+aftercut*(n-1):aftercut*n,i)= bufi(start+beforecut*(n-1):over+beforecut*(n-1),i);
        end
    end
    plot(cut(:,i) - maxsample*(i-1)*1,'LineWidth',.6)
    hold on
end
yticks([-(maxsample*(channel_num-1)):maxsample:0])
for i=1:channel_num
    ticklabels{i} = strcat('ch',num2str(channel_num-i+1));
end
yticklabels(ticklabels)
xlim([-inf,inf])
ylim([-(maxsample*channel_num),maxsample])
ylabel('Channels')
xlabel('Samples')
handles.cut = cut;
guidata(hObject, handles)
now = cut;
handles.now = now;
guidata(hObject,handles)

function pushbutton9_Callback(hObject, eventdata, handles)
% clean axes2
try
    delete(allchild(handles.axes2));
end

for i = 1:4
    olddata = get(handles.uitable1,'Data');%获得目前表格里的内容
    row = size(olddata,1);
    tempdata = olddata(1:row-1,:);
    set(handles.uitable1,'Data',tempdata);
end

channel_num = handles.channel_num;
gesture_num = handles.gesture_num;
repeat = handles.repeat;
cut = handles.cut;
aftercut = handles.aftercut;
maxsample = handles.maxsample;

% plot normolizetion
axes(handles.axes2);  
for i=1:channel_num
    for k=1:gesture_num
        for j = 1:repeat
            n = j+(k-1)*repeat;
            c = cut(1+aftercut*(n-1):aftercut*n,i);
            temp =[];
            demesion = size(c);
            mean_nor = mean(c);
            std_nor = std(c);
            for p = 1:numel(c)
                temp = [temp (c(p)-mean_nor)/std_nor];
            end
            nor(1+aftercut*(n-1):aftercut*n,i) = reshape(temp,demesion);
        end
    end
    plot(nor(:,i) - maxsample*(i-1)*1,'LineWidth',.6)
    hold on
end
handles.nor = nor;
guidata(hObject, handles)
now = nor;
handles.now = now;
guidata(hObject, handles)
yticks([-(maxsample*(channel_num-1)):maxsample:0])
for i=1:channel_num
    ticklabels{i} = strcat('ch',num2str(channel_num-i+1));
end
yticklabels(ticklabels)
xlim([-inf,inf])
ylim([-(maxsample*channel_num),maxsample])
ylabel('Channels')
xlabel('Samples')

function checkbox1_Callback(hObject, eventdata, handles)

function checkbox2_Callback(hObject, eventdata, handles)

function checkbox3_Callback(hObject, eventdata, handles)

function checkbox4_Callback(hObject, eventdata, handles)

function checkbox5_Callback(hObject, eventdata, handles)

function checkbox6_Callback(hObject, eventdata, handles)

function checkbox7_Callback(hObject, eventdata, handles)

function checkbox8_Callback(hObject, eventdata, handles)

function checkbox9_Callback(hObject, eventdata, handles)

function checkbox10_Callback(hObject, eventdata, handles)

function checkbox11_Callback(hObject, eventdata, handles)

function checkbox12_Callback(hObject, eventdata, handles)

function radiobutton1_Callback(hObject, eventdata, handles)

function radiobutton2_Callback(hObject, eventdata, handles)

function radiobutton3_Callback(hObject, eventdata, handles)

function radiobutton4_Callback(hObject, eventdata, handles)

function radiobutton5_Callback(hObject, eventdata, handles)

function radiobutton6_Callback(hObject, eventdata, handles)

function radiobutton7_Callback(hObject, eventdata, handles)

function pushbutton5_Callback(hObject, eventdata, handles)
axes(handles.axes3); %以下作图在axes3
nor = handles.nor;
channel_num = handles.channel_num;
aftercut = handles.aftercut;
fs = handles.fs;
gesture_num = handles.gesture_num;
repeat = handles.repeat;

v1 = get(handles.checkbox1,'Value');
v2 = get(handles.checkbox2,'Value');
v3 = get(handles.checkbox3,'Value');
v4 = get(handles.checkbox4,'Value');
v5 = get(handles.checkbox5,'Value');
v6 = get(handles.checkbox6,'Value');
v7 = get(handles.checkbox7,'Value');
v8 = get(handles.checkbox8,'Value');
v9 = get(handles.checkbox9,'Value');
v10 = get(handles.checkbox10,'Value');
v11 = get(handles.checkbox11,'Value');
v12 = get(handles.checkbox12,'Value');
c1 = get(handles.radiobutton1,'Value');
c2 = get(handles.radiobutton2,'Value');
c3 = get(handles.radiobutton3,'Value');
c4 = get(handles.radiobutton4,'Value');
c5 = get(handles.radiobutton5,'Value');

% Feature Extraction
check = 0;  % number of chosen checkboxes
su = 0;  % number of the features

% Sliding window
f_len = floor(aftercut*0.6);
f_num = 10;
f_step = floor((aftercut-f_len)/(f_num-1));

each_feature_num = gesture_num*repeat*f_num;

% 1-RMS
if(v1 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(1);
    else
        fet = strcat('+',num2str(1));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    feature(each_feature_num*su+fm,i) = sqrt(mean(data.^2));
                end
            end
        end
    end
    su = su+1;
end

% 2-WL
if(v2 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(2);
    else
        fet = strcat('+',num2str(2));
        fe = strcat(fe,fet);
        disp(fe);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    feature(each_feature_num*su+fm,i) = sum(abs(diff(data)))/length(data);
                end
            end
        end
    end
    su = su+1;
end

% 3-MAV
if(v3 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(3);
    else
        fet = strcat('+',num2str(3));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    feature(each_feature_num*su+fm,i) = mean(abs(data));
                end
            end
        end
    end
    su = su+1;
end

% 4-VAR
if(v4 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(4);
    else
        fet = strcat('+',num2str(4));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    feature(each_feature_num*su+fm,i) = var(data);
                end
            end
        end
    end
    su = su+1;
end

% 5-ZC
if(v5 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(5);
    else
        fet = strcat('+',num2str(5));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    deadzone = 10e-7;
                    data_size = length(data);
                    feature(each_feature_num*su+fm,i) = 0;
                    if data_size == 0
                        feature(each_feature_num*su+fm,i) = 0;
                    else
                        for q = 2:data_size
                            difference = data(q)-data(q-1);
                            multy = data(q)*data(q-1);
                            if abs(difference)>deadzone && multy<0
                                feature(each_feature_num*su+fm,i) = feature(each_feature_num*su+fm,i)+1;
                            end
                        end
                        feature(each_feature_num*su+fm,i) = feature(each_feature_num*su+fm,i)/data_size;
                    end
                end
            end
        end
    end
    su = su+1;
end

% 6-SSC
if(v6 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(6);
    else
        fet = strcat('+',num2str(6));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    deadzone = 10e-7;
                    data_size = length(data);
                    feature(each_feature_num*su+fm,i) = 0;
                    if data_size == 0
                        feature(each_feature_num*su+fm,i) = 0;
                    else
                        for q = 3:data_size
                            difference1 = data(q-1)-data(q-2);
                            difference2 = data(q-1)-data(q);
                            sign = difference1*difference2;
                            if sign>0
                                if abs(difference1)>deadzone || abs(difference2)>deadzone
                                    feature(each_feature_num*su+fm,i) = feature(each_feature_num*su+fm,i)+1;
                                end
                            end
                        end
                        feature(each_feature_num*su+fm,i) = feature(each_feature_num*su+fm,i)/data_size;
                    end
                end
            end
        end
    end
    su = su+1;
end

% 7-MNF
if(v7 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(7);
    else
        fet = strcat('+',num2str(7));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    feature(each_feature_num*su+fm,i) = meanfreq(data,fs);
                end
            end
        end
    end
    su = su+1;
end

% 8-MDF
if(v8 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(8);
    else
        fet = strcat('+',num2str(8));
        fe = strcat(fe,fet);
    end
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    feature(each_feature_num*su+fm,i) = medfreq(data,fs);
                end
            end
        end
    end
    su = su+1;
end

% 9-Haar(db1)
if(v9 == 1)
    check = check+1;
    if(check == 1)
        fe = num2str(9);
    else
        fet = strcat('+',num2str(9));
        fe = strcat(fe,fet);
    end   
    
    wavename = 'haar';
    N = 4;  % order    
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    [C,L] = wavedec(data, N, wavename);
                    d4 = wrcoef('d',C,L,wavename,4);  
                    a4 = wrcoef('a',C,L,wavename,4);
                    % get coef's RMS as features
                    feature(each_feature_num*su+fm,i) = norm(d4);
                    feature(each_feature_num*(su+1)+fm,i) = norm(a4);
                end
            end
        end
    end
    su = su+2;
end

% 10-CoifN
if(v10 == 1)
    check = check+1;
    n_coif = get(handles.edit12,'string');
    if(check == 1)
        fe = strcat(num2str(10),'(');
        fe = strcat(fe,n_coif);
        fe = strcat(fe,')');
    else
        fet = strcat('+',num2str(10));
        fe = strcat(fe,fet);
        fe = strcat(fe,'(');
        fe = strcat(fe,n_coif);
        fe = strcat(fe,')');
    end

    wavename = strcat('coif',n_coif);
    N = 4;
    
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    [C,L] = wavedec(data, N, wavename);
                    d4 = wrcoef('d',C,L,wavename,4);
                    a4 = wrcoef('a',C,L,wavename,4);
                    feature(each_feature_num*su+fm,i) = norm(d4);
                    feature(each_feature_num*(su+1)+fm,i) = norm(a4);
                end 
            end
        end
    end
    su = su+2;
end

% 11-DbN
if(v11 == 1)
    check = check+1;    
    n_db = get(handles.edit8,'string');    
    if(check == 1)
        fe = strcat(num2str(11),'(');
        fe = strcat(fe,n_db);
        fe = strcat(fe,')');
    else
        fet = strcat('+',num2str(11));
        fe = strcat(fe,fet);
        fe = strcat(fe,'(');
        fe = strcat(fe,n_db);
        fe = strcat(fe,')');
    end
    
    wavename = strcat('db',n_db);
    N = 4;
    
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    [C,L] = wavedec(data, N, wavename);
                    d4 = wrcoef('d',C,L,wavename,4);
                    a4 = wrcoef('a',C,L,wavename,4);
                    feature(each_feature_num*su+fm,i) = norm(d4);
                    feature(each_feature_num*(su+1)+fm,i) = norm(a4);
                end 
            end
        end
    end
    su = su+2;
end

% 12-SymN
if(v12 == 1)
    check = check+1;
    n_sym = get(handles.edit9,'string');%从edit9获取symN的N    
    if(check == 1)
        fe = strcat(num2str(12),'(');
        fe = strcat(fe,n_sym);
        fe = strcat(fe,')');
    else
        fet = strcat('+',num2str(12));
        fe = strcat(fe,fet);
        fe = strcat(fe,'(');
        fe = strcat(fe,n_sym);
        fe = strcat(fe,')');
    end
    
    wavename = strcat('sym',n_sym);
    N = 4;
    
    for i=1:channel_num
        for k=1:gesture_num
            for j = 1:repeat
                n = j+(k-1)*repeat;
                for f = 1:f_num
                    fm = (n-1)*f_num+f;
                    f_start = 1+aftercut*(n-1)+f_step*(f-1);
                    f_end = f_start+f_len;
                    data = nor(f_start:f_end,i);
                    [C,L] = wavedec(data, N, wavename);
                    d4 = wrcoef('d',C,L,wavename,4);
                    a4 = wrcoef('a',C,L,wavename,4);
                    feature(each_feature_num*su+fm,i) = norm(d4);
                    feature(each_feature_num*(su+1)+fm,i) = norm(a4);
                end 
            end
        end
    end
    su = su+2;
end

% Classification  
eachnum = [1:repeat*f_num];
t1 = round(repeat*f_num*0.75);  % 75 percent for train
t2 = repeat*f_num-t1;
random_num = eachnum(randperm(numel(eachnum),t1));
train_num = sort(random_num);
isb = ismember(eachnum,train_num);
test_num = eachnum(~isb);

% trainset and testset
for fe_num = 1:su
    x_all(1:repeat*gesture_num*f_num,1+channel_num*(fe_num-1):channel_num*fe_num) = feature(1+repeat*gesture_num*f_num*(fe_num-1):f_num*repeat*gesture_num*fe_num,1:channel_num);
end
handles.x_all = x_all;
guidata(hObject, handles)
for i = 1:gesture_num
    y_train(1+t1*(i-1):t1*i,1) = i;
    y_test(1+t2*(i-1):t2*i,1) = i;
end
for i = 1:gesture_num
    for j = 1:t1
        x_train(j+t1*(i-1),:) = x_all(train_num(j)+f_num*repeat*(i-1),:);
    end
    for p = 1:t2
        x_test(p+t2*(i-1),:) = x_all(test_num(p)+f_num*repeat*(i-1),:);
    end
end

% 1-SVM
% download libsvm toolbox
if(c1 == 1)
    cf = num2str(1);
    [x_train,pstrain] = mapminmax(x_train');
    pstrain.ymin = 0;
    pstrain.ymax = 1;
    [x_train,pstrain] = mapminmax(x_train,pstrain);  % normalization [0,1]
    [x_test,pstest] = mapminmax(x_test');
    pstest.ymin = 0;
    pstest.ymax = 1;
    [x_test,pstest] = mapminmax(x_test,pstest);
    x_train = x_train';
    x_test = x_test';
    
    % find best c and g
    % c&g approximately range 2^(-10),2^(-9),...,2^(10)
    [bestacc,bestc,bestg] = SVMcgForClass(y_train,x_train,-10,10,-10,10);
    % c precise 2^(-2),2^(-1.5),...,2^(4) and g 2^(-4),2^(-3.5),...,2^(4)
    [bestacc,bestc,bestg] = SVMcgForClass(y_train,x_train,-2,4,-4,4,3,0.5,0.5,0.9);
    % train for model
    cmd = ['-c ',num2str(bestc),' -g ',num2str(bestg)];
    model=libsvmtrain(y_train,x_train,cmd);   
    % predict
    [predict_label, accuracy3, dec_values]=libsvmpredict(y_test,x_test,model);
    % accuracy3(3*1):accuracy、rmserror、squared multiple correlation coefficient
    accuracy = accuracy3(1,1);
end

% 2-KNN
if(c2 == 1)
    cf = num2str(2);
    k_knn = str2double(get(handles.edit10,'string')); 
    cf = strcat(cf,'(');
    cf = strcat(cf,get(handles.edit10,'string'));
    cf = strcat(cf,')');
    
    mdl = fitcknn(x_train,y_train,'NumNeighbors',k_knn)
    class = predict(mdl,x_test);
    count = 0;
    for n = 1:gesture_num*t2
        if(class(n) == y_test(n))
            count = count+1;
        end
    end
    accuracy = count/(gesture_num*t2);
    accuracy = accuracy*100;
end   

% 3-LDA
if(c3 == 1)
    cf = num2str(3);
    %w_i为第i类样本
    %W为权值向量
    %第一步：计算样本均值向量
    for i = 1:gesture_num
        w_i = x_train(1+t1*(i-1):t1*i,:);
        M_1(i,:) = mean(w_i);%求每列均值,M是行向量,列数为su*channel_num
    end
    M = mean(x_train);
    % Sw
    Sw = zeros(su*channel_num);
    for g = 1:gesture_num
        S_g = zeros(su*channel_num);
        w_g = x_train(1+t1*(g-1):t1*g,:);
        for i = 1:t1
            S_g = S_g + (w_g(i,:)-M_1(g,:))'*(w_g(i,:)-M_1(g,:));
        end     
        Sw = Sw + S_g;
    end
    % Sb
    Sb = zeros(su*channel_num);
    for g = 1:gesture_num
        Si = M - M_1(g,:);
        Sb = Sb + Si'*Si;
    end
    % find max feature vector W
    [V,L]=eig(inv(Sw)*Sb);
    [a,b]=max(max(L));
    W = V(:,b);   
    for g = 1:gesture_num
        w_g = x_train(1+t1*(g-1):t1*g,:);
        for i = 1:t1
            y(i+t1*(g-1)) = W'* w_g(i,:)';
        end
        m(g,:) = mean(y(1+t1*(g-1):t1*g));
    end
    % predict
    correct = 0;
    for j = 1:t2*gesture_num
        test_y=W'*x_test(j,:)';
        min_index=0;
        min_dis=2;
        for i=1:gesture_num
            if(abs(test_y-m(i,:)))<min_dis
                min_index=i;
                min_dis=abs(test_y-m(i,:));
            end
        end
        if min_index == y_test(j,1)
            correct = correct+1;
        end
    end
    accuracy = (correct/(t2*gesture_num))*100;
end

% 4-ANN
if(c4 == 1)
    cf = num2str(4);
    x_train = x_train';
    s = length(y_train) ;
    output = zeros(s,gesture_num) ;
    for i = 1:s
        output(i,y_train(i)) = 1;
    end
    output_train = output';
    
    [input_train_m,ps] = mapminmax(x_train,0,1);  
    num_layer = round(sqrtm(su*channel_num+gesture_num)+10);
    net = newff(minmax(input_train_m),[100,num_layer,gesture_num],{'logsig','logsig','logsig'},'traingdx');
    net.trainParam.min_grad=1e-25;     
    net.trainparam.epochs=10000;        
    net.trainparam.goal=0.0000001;      
    net.trainParam.lr=0.01;             
    net = train(net,input_train_m,output_train);
    input_test=mapminmax('apply',x_test',ps);
    % predict
    testoutput=sim(net,input_test);
    [s1,s2] = size(testoutput);
    right_number = 0 ;
    for i = 1 : s2
        [c,Index]= max(testoutput(:,i));
        if(Index == y_test(i))
            right_number = right_number + 1 ;
        end
    end
    accuracy = (right_number/s2)*100;
end

% 5-Random Forest
if(c5 == 1) 
    cf = num2str(5);
    nTree=10;
    B=TreeBagger(nTree,x_train,y_train,'Method', 'classification');
    predictl=predict(B,x_test);
    class=str2double(predictl);
    count = 0;
    for n = 1:gesture_num*t2
        if(class(n) == y_test(n))
            count = count+1;
        end
    end
    accuracy = (count/(gesture_num*t2))*100;
end
      
olddata = get(handles.uitable1,'Data');%获得目前表格里的内容
accuracy = roundn(accuracy,-3);
newdata = {fe,cf,num2str(accuracy)};
tempdata = [olddata;newdata];
set(handles.uitable1,'Data',tempdata);

function pushbutton6_Callback(hObject, eventdata, handles)
fe = 'blank';
nor = handles.nor;
aftercut = handles.aftercut;
channel_num = handles.channel_num;
gesture_num = handles.gesture_num;
repeat = handles.repeat;
f_len = floor(aftercut*0.6);
f_num = 10;
f_step = floor((aftercut-f_len)/(f_num-1));

for k=1:gesture_num
    for j = 1:repeat
        n = j+(k-1)*repeat;
        for f = 1:f_num
            fm = (n-1)*f_num+f;
            f_start = aftercut*(n-1)+1+f_step*(f-1);
            f_end = aftercut*(n-1)+f_len+f_step*(f-1);
            x_all(1+f_len*(fm-1):f_len*fm,:) = nor(f_start:f_end,:);
            y_all(fm,1) = k;
        end
    end
end

% trainset and testset
eachnum = [1:repeat*f_num];
t1 = round(repeat*f_num*0.75);
t2 = repeat*f_num-t1;
random_num = eachnum(randperm(numel(eachnum),t1));
train_num = sort(random_num);
isb = ismember(eachnum,train_num);
test_num = eachnum(~isb);
for i = 1:gesture_num
    for j = 1:t1
        n1 = j+t1*(i-1);
        m1 = train_num(j)+f_num*repeat*(i-1);
        x_train(1+f_len*(n1-1):f_len*n1,:) = x_all(1+f_len*(m1-1):f_len*m1,:);
        y_train(n1,1) = i;
    end
    for p = 1:t2
        n2 = p+t2*(i-1);
        m2 = test_num(p)+f_num*repeat*(i-1);
        x_test(1+f_len*(n2-1):f_len*n2,:) = x_all(1+f_len*(m2-1):f_len*m2,:);
        y_test(n2,1) = i;
    end
end

c6 = get(handles.radiobutton6,'Value');
c7 = get(handles.radiobutton7,'Value');

% 6-CNN
if(c6 ==  1)
    cf = num2str(6);
    savepath = pwd;
    file_path_name = strcat(savepath,'\');
    file_path_name = strcat(file_path_name,'signal_pic');
    path_x = file_path_name;
    file_path = strcat(file_path_name,'\');    
    for i=1:gesture_num
        file_path_name = strcat(file_path,num2str(i));
        
        if exist(file_path_name)==0 
            mkdir(file_path_name);
        else 
            rmdir(file_path_name, 's'); %该文件夹中有没有文件均可
            mkdir(file_path_name);
        end
        
        for j = 1:repeat*f_num
            p=repeat*f_num*(i-1)+j;
            data = x_all(1+f_len*(p-1):f_len*p,:)';
            path = strcat(file_path_name,'\');           
            path = strcat(path,num2str(j));
            path = strcat(path,'.jpg');           
            imwrite(data,path);
        end
    end
    
    imds = imageDatastore(path_x,'IncludeSubfolders',true,'FileExtensions','.jpg',...
        'LabelSource','foldernames');
    labelCount = countEachLabel(imds);
    img = readimage(imds,1);
    size(img);
    numTrainFiles = t1;
    [imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');
    
    layers = [
        imageInputLayer([channel_num f_len 1])
        
        convolution2dLayer(5,20,'Padding','same')
        batchNormalizationLayer
        reluLayer
        
        maxPooling2dLayer(2,'Stride',2)
        
        convolution2dLayer(5,16,'Padding','same')
        batchNormalizationLayer
        reluLayer
        
        maxPooling2dLayer(2,'Stride',2)
        
        convolution2dLayer(3,32,'Padding','same')
        batchNormalizationLayer
        reluLayer
        
        fullyConnectedLayer(gesture_num)
        softmaxLayer
        classificationLayer];
    
    options = trainingOptions('sgdm', ...
        'InitialLearnRate',0.001, ...
        'MaxEpochs',100, ...
        'Shuffle','every-epoch', ...
        'ValidationData',imdsValidation, ...
        'ValidationFrequency',30, ...
        'Verbose',false, ...
        'Plots','training-progress');
    net = trainNetwork(imdsTrain,layers,options);
    YPred = classify(net,imdsValidation);
    YValidation = imdsValidation.Labels;
    
    accuracy = (sum(YPred == YValidation)/numel(YValidation))*100;
end

% 7-RNN
if(c7 ==  1)
    cf = num2str(7);
    XTrainData = x_train';
    XTrainLabel = y_train;
    %XTrain
    for i = 1:t1*gesture_num
        XTrain{i,1} = XTrainData(:,1+f_len*(i-1):f_len*i);
    end
    %YTrain
    TrainstrLabel = num2str(XTrainLabel);  % num to str
    
    for i = 1:t1*gesture_num  % str matrix to cell
        TraincellLabel{i,1} = TrainstrLabel(i,:);
        
    end
    YTrain = categorical(TraincellLabel);  % cell to categorical
    % testset part
    XTestData = x_test';
    XTestLabel = y_test;
    for i = 1:t2*gesture_num
        XTest{i,1} = XTestData(:,1+f_len*(i-1):f_len*i);
    end
    TeststrLable = num2str(XTestLabel);  % num to str
    for i = 1:t2*gesture_num
        TestcellLable{i,1} = TeststrLable(i,:);  % str matrix to cell
    end
    YTest = categorical(TestcellLable);  % cell to categorical
    %layers
    inputSize = channel_num;  %将输入大小指定为序列大小 10（输入数据的维度,指同一时间下的数据维度）
    numHiddenUnits = 100;  %指定具有 100 个隐含单元的双向 LSTM 层，并输出序列的最后一个元素
    numClasses = gesture_num;  %指定12个类，包含大小为 1 的全连接层，后跟 softmax 层和分类层
    
    layers = [ ...
        sequenceInputLayer(inputSize)
        bilstmLayer(numHiddenUnits,'OutputMode','last')
        fullyConnectedLayer(numClasses)
        softmaxLayer
        classificationLayer];
    %sequenceInputLayer(inputSize)：序列输入层，指定输入维度
    %bilstmLayer(numHiddenUnits,'OutputMode','last')：双向LSTM层，指定隐藏节点，输出模式为‘last’即输出最后一个分类值
    %fullyConnectedLayer(numClasses)：全连接层，指定输出类别的个数
    %softmaxLayer：这层是输出各类别分类的的概率
    %classificationLayer：分类层，输出最后的分类结果，类似于概率竞争投票
    
    %options
    maxEpochs = 100;
    miniBatchSize = 50;
    
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'GradientThreshold',1, ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'SequenceLength','longest', ...
        'Shuffle','never', ...
        'Verbose',0, ...
        'Plots','training-progress');
    %求解器为 'adam'
    %梯度阈值为 1
    %最大轮数为 100
    %50作为小批量数
    %填充数据以使长度与最长序列相同，序列长度指定为 'longest'。
    %数据保持按序列长度排序的状态，不打乱数据。
    %'ExecutionEnvironment' 指定为 'cpu'，设定为'auto'表明使用GPU。
    
    %训练网络
    net = trainNetwork(XTrain,YTrain,layers,options);
    
    %利用网络来分类
    miniBatchSize = 30;
    YPred = classify(net,XTest, ...
        'SequenceLength','longest','MiniBatchSize',miniBatchSize);
    %计算分类准确度
    accuracy = (sum(YPred == YTest)./numel(YTest))*100;
end

%表格部分
olddata = get(handles.uitable1,'Data');
accuracy = roundn(accuracy,-3);
newdata = {fe,cf,num2str(accuracy)};
tempdata = [olddata;newdata];
set(handles.uitable1,'Data',tempdata);

function pushbutton7_Callback(hObject, eventdata, handles)
derow = str2double(get(handles.edit11,'string'));
olddata = get(handles.uitable1,'Data');
olddata(derow,:) = [];
set(handles.uitable1,'Data',olddata);

function pushbutton12_Callback(hObject, eventdata, handles)
% draw the bar
axes(handles.axes3);
olddata = get(handles.uitable1,'Data');
nr = size(olddata);
x = [1:1:nr];
for n_r = 1:nr
    y(n_r) = str2num(olddata{n_r,3});
end
b = bar(x,y,0.5);
set(b,'FaceColor',[1,0.9,0.3]);
for n_r = 1:nr
    text(x(n_r),y(n_r),num2str(y(n_r)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8);
end
xlabel('Index');
ylabel('Accuracy');
set(gca,'ylim',[0,100],'ytick',[0:10:100]);

function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit10_Callback(hObject, eventdata, handles)

function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)

function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles)

function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_SizeChangedFcn(hObject, eventdata, handles)

function pushbutton10_Callback(hObject, eventdata, handles)
now = handles.now;
your_data = now;
[filename filepath]=uiputfile('*.mat');
savefilename = [filepath,filename];
save(savefilename,'your_data');

function pushbutton11_Callback(hObject, eventdata, handles)
x_all = handles.x_all;
your_feature = x_all;
[filename filepath]=uiputfile('*.mat');
savefilename = [filepath,filename];
save(savefilename,'your_feature');
