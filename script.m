function app
%% Info %%

% Author: Christian Menges, 2021

% This script analyzes data produced by Equa IDA ICE or other tools with
% similar output formats. It produces distribution plots, cumulative plots,
% statistical plots, etc., and calculates degree/ppm hours. In general, the
% script is independent of the type of input data, but the main
% application areas are the evaluation of temperature and CO2 profiles.

%% Usage %%

% Adapt the settings below to the properties of your simulation.
% Afterwards, execute the script by pressing 'F5' or by clicking the button
% 'Run' in the Matlab programming environment. Depending on the settings
% you will be prompted to select the input file (.xlxs (Excel file)), the
% sheet name inside of the input file as well as the data column, which
% should be analyzed. If the setting 'useIdealTemperature' is set to
% 'true', an additional file specifying the weather conditions must be
% selected. Weather conditions are assumed to be stored as .prn file there
% the first column represents the passed time in hours relative to the
% beginning of the simulation (usually the beginning of the year) and the
% second column shows the environment temperature in degree Celsius. Once
% all files are selected, several windows will show up presenting you the
% analysis results. Additionally, most of the results are written into an
% Excel file, which can be used for further data exploration.

% Note: As of May 2021, IDA ICE doesn't correctly handle daylight saving
% time. This leads to an offset in the data, which can be visualized by
% using the carpet plot feature of IDA ICE. This script does account for
% this bug and corrects the wrong data.
 
% For questions or bug reports, please contact: christian.menges@tum.de

%% Configuration %%
simulationStart = '2009-01-01 00:00:00'; % start of the simulation

% Select time frame for analysis
% (endDate = startDate uses all data after startdate)
startDate = simulationStart;
endDate = simulationStart;

% Sheet name in the input file, if empty a selection dialog is shown
% example: 
% "Indoor air quality measures_ Ac" or "Air and operative temperatures"
sheetName = ""; 

% Column name in the data sheet in input file, if empty a selection dialog
% is shown
columnName = ""; % common values: "CO2_Ppm_vol_";

% Use a threshold value depending on the environment temperature
useIdealTemperature = true; % false
% .prn file for climate conditions
climateFile = ""

% Type and unit of analyzation (used for data labels)
type = 'degree'; % 'CO2'
unit = '°C'; % 'ppm'

% If relative analysis is wanted set to CO2 level / temperature of the
% outside air, else set to 0
relativeLevel = 0; 

% Schedule representing time of interest, often usage time
% (if complete days should be taken into account use [0 24])
% Note: Currently, only full hours are supported!
MoFr = [7 18];
Saturday = [];
Sunday = [];

% Holidays without usage
% Input format: Matrix, first column represents the start and the
% second column the end of each holiday. 
% example:
% holidays = [
%     "2009-12-24 00:00:00" "2009-12-31 23:59:00";
%     ...
% ]
holidays = []

% holidays = [
%     "2009-12-24 00:00:00" "2009-12-31 23:59:00";
%     "2009-11-18 00:00:00" "2009-11-18 23:59:00";
%     "2009-11-02 00:00:00" "2009-11-03 23:59:00";
%     "2009-08-03 00:00:00" "2009-09-14 23:59:00";
%     "2009-06-01 00:00:00" "2009-06-12 23:59:00";
%     "2009-05-21 00:00:00" "2009-05-21 23:59:00";
%     "2009-05-01 00:00:00" "2009-05-01 23:59:00";
%     "2009-04-06 00:00:00" "2009-04-17 23:59:00";
%     "2009-02-23 00:00:00" "2009-02-27 23:59:00";
%     "2009-01-01 00:00:00" "2009-01-06 23:59:00"
%     ];

% Ranges in bar plot
%edges = [0,1000,100000];
%edges = [0,850,1000,2000,3000,100000];
edges = [0, 400:100:1300, 10000];

% Threshold value to calculate CO2/degree hours (ignored if ideal
% temperature is used)
threshold = 1000;

% Name of the exported Excel file
exportFilename = 'results.xlsx';

interpolationMethod = 'spline'; % 'spline', 'linear'
sampleInterval = 6; % Minutes

%% Program %%
DateInputFormat = 'yyyy-MM-dd HH:mm:ss';
simulationStartInternal = convertDate(simulationStart, DateInputFormat);
startDateInternal = convertDate(startDate,DateInputFormat);
endDateInternal = convertDate(endDate,DateInputFormat);

if(endDateInternal < startDateInternal)
    error('Error: endDate < startDate');
end

if(startDateInternal < simulationStartInternal)
    warning('startDate < simulationStartInternal');
end

f = figure('Visible','off','Position',[500,500,450,285]);
f.Name = strcat(type, ' Evaluation');
movegui(f,'center');

% Open input file
filePath = selectFile('*.xlsx', 'Select input data file (.xlsx)');

if sheetName == ""
    sheetName = selectOption(sheetnames(filePath), 'Select data sheet');
end

T = readtable(filePath, 'Sheet', sheetName);

if columnName == ""
    columnName = char(selectOption(T.Properties.VariableNames, 'Select data column'));
end

if ~any(strcmp('Time', T.Properties.VariableNames))
   if any(strcmp('Zeit', T.Properties.VariableNames))
       T.Properties.VariableNames{'Zeit'} = 'Time';
   else
       error('Error: Input data sheet does not contain column named "Time" or "Zeit"');
   end
end

querypoints = 0:sampleInterval/60:max(T.Time);

outsideTemperatur = zeros(1, length(querypoints)) + threshold;
if useIdealTemperature
    if climateFile == ""
        climateFile = selectFile('*.prn', 'Select climate file (.prn)');
    end
    % Read outside temperature
    outsideTemperatur = readmatrix(climateFile, 'FileType', 'text');
    outsideTemperatur = outsideTemperatur(:,[1,2]);
end

% Convert input data
T = sortrows(T);
[~, ia, ~] = unique(T.Time, 'rows', 'sorted');
values = interp1(T(ia,:).Time, table2array(T(ia, columnName)), querypoints, interpolationMethod);
dates = simulationStartInternal + datenum(hours(querypoints));
datesTime = datetime(dates.', 'ConvertFrom', 'datenum', 'TimeZone', 'Europe/Zurich');
[weekdayNum, weekdayStr] = weekday(datesTime);

weekdayNum = weekdayNum.';
values = values.';
values = values - relativeLevel;
if ~useIdealTemperature && strcmp(type, 'CO2') == 1
    values(values < 0) = 0;
end

if useIdealTemperature
    outsideTemperatur = interp1(outsideTemperatur(:,1), outsideTemperatur(:,2), querypoints, interpolationMethod);
end

% data table format
% | weekdayNum | weekdayStr | daylight saving time (dst) | dates | year |
% month | day | hour | minute | sec | threshold |
data = [table(weekdayNum.',string(weekdayStr), isdst(datesTime), dates.', 'VariableNames', {'weekdayNum', 'weekdayStr', 'dst','dates'}), ...
    array2table(datevec(datenum(datesTime)+datenum(hours(isdst(datesTime)))), 'VariableNames', {'Year', 'Month', 'Day', 'Hour', 'Minute','Sec'}), ...
    table(values),table(outsideTemperatur.', 'VariableNames', {'threshold'})];

data = filterData(data, simulationStart, startDate, endDate, MoFr, Saturday, Sunday, holidays, DateInputFormat);

% Analyze data and plot results
hist = histcounts(data.values, edges);
hist = hist / (60 / sampleInterval);

xlabels = createXLabels(edges);

plotBarPlot(hist, xlabels, unit, type);

if useIdealTemperature
    [RhMonthsDIN20, RunithMonthsDIN20, Runith2MonthsDIN20] = calculateRevUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureDIN(e,20));
    
    [hMonthsDIN20, unithMonthsDIN20, unith2MonthsDIN20] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureDIN(e,20));
    [hMonthsDIN22, unithMonthsDIN22, unith2MonthsDIN22] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureDIN(e,22));
    [hMonthsDIN24, unithMonthsDIN24, unith2MonthsDIN24] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureDIN(e,24));
    
    % Plot results
    figure, plot(1:12, unithMonthsDIN20, 1:12, unithMonthsDIN22,'--', 1:12, unithMonthsDIN24), 
    title('DIN EN 15251:2012-12'), xlabel("Month"), ylabel(strcat("Delta ", unit, " hours"));
    figure, plot(1:12, unith2MonthsDIN20, 1:12, unith2MonthsDIN22,'--', 1:12, unith2MonthsDIN24), 
    title('DIN EN 15251:2012-12 squared'),xlabel("Month"), ylabel(strcat("Delta² ", unit, " hours"));
    
    [RhMonthsASHRAE143, RunithMonthsASHRAE143, Runith2MonthsASHRAE143] = calculateRevUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureASHRAE(e,14.3));
    
    [hMonthsASHRAE143, unithMonthsASHRAE143, unith2MonthsASHRAE143] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureASHRAE(e,14.3));
    [hMonthsASHRAE178, unithMonthsASHRAE178, unith2MonthsASHRAE178] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureASHRAE(e,17.8));
    [hMonthsASHRAE213, unithMonthsASHRAE213, unith2MonthsASHRAE213] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, @(e) idealTemperatureASHRAE(e,21.3));
    
    % Plot results
    figure, plot(1:12, unithMonthsASHRAE143, 1:12, unithMonthsASHRAE178,'--', 1:12, unithMonthsASHRAE213),
    title('ASHRAE 55-2020'), xlabel("Month"), ylabel(strcat("Delta ", unit, " hours"));
    figure, plot(1:12, unith2MonthsASHRAE143, 1:12, unith2MonthsASHRAE178,'--', 1:12, unith2MonthsASHRAE213), 
    title('ASHRAE 55-2020 squared'), xlabel("Month"), ylabel(strcat("Delta² ", unit, " hours"));
else
    [hMonths, unithMonths, unith2Months] = calculateUnitHours(data, sampleInterval, unit);
    
    % Plot results
    figure, plot(1:12, unithMonths), title('ASHRAE 55-2020'), 
    xlabel("Month"), ylabel(strcat("Delta ", unit, " hours"));
    figure, plot(1:12, unith2Months), title('ASHRAE 55-2020'), 
    xlabel("Month"), ylabel(strcat("Delta² ", unit, " hours"));
end

usageTimeMonths = zeros(12,1);
for m = 1:12
    usageTimeMonths(m) = height(data.values(data.Month == m)) / (60 / sampleInterval);
end

if useIdealTemperature == false
    T2 = data.values(data.values > data.threshold) - data.threshold(data.values > data.threshold);
    T2 = sortrows(T2);
    S = cumsum(T2/ (60 / sampleInterval));
    figure, plot(T2,S), xlabel(strcat("Delta ", unit)), ylabel(strcat("Delta ", unit, " hours"));
end

figure, subplot('position', [0.05 0.1 0.7 0.7]),boxplot(data.values, month(datetime(1, data.Month, 1), 'name')),
xlabel("Months"), ylabel(unit);
subplot('position', [0.85 0.1 0.1 0.7]), boxplot(data.values),
xlabel(strcat(datestr(startDateInternal, 'dd.mm.yyyy'), " - ", datestr(endDateInternal, 'dd.mm.yyyy'))), ylabel(unit);

f.Visible = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excel export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfile(exportFilename)
   answer = questdlg(strcat(exportFilename, " already exists. Overwrite file?"), 'Warning', 'Yes', 'Cancel', 'Yes')
   if strcmp(answer, 'Cancel')
       return;
   end
   % Remove old file, otherwise data can mix
   delete(exportFilename);
end

T = table(xlabels.', hist.', 'VariableNames', {'Range', 'Hours'});
writetable(T, exportFilename, 'Sheet', 'Hours');
[meanValues, minValues, maxValues, valueCount] = grpstats(data.values, data.Month, {'mean','min','max','numel'}, 'Alpha', 0.25);
% Pass function to grpstats as described here:
% https://de.mathworks.com/matlabcentral/answers/156729-how-to-use-prctile-within-grpstats#answer_153537
percentiles = grpstats(data.values, data.Month, @(y) prctile(y, [5 25 50 75 95]'));
months = unique(data.Month);
varNames = {'percentile5','percentile25','percentile50','percentile75','percentile95'};
T = [table(months, meanValues, minValues, maxValues) array2table(percentiles, 'VariableNames', varNames) table(valueCount)];
writetable(T, exportFilename, 'Sheet', 'Statistics');

if useIdealTemperature == false
T = table(T2, S, 'VariableNames', {'Diff', 'Sum'});
writetable(T, exportFilename, 'Sheet', 'Sum');
end

m = 1:12;
if useIdealTemperature
    [minValuesT, maxValuesT] = grpstats(data.threshold, data.Month, {'min','max'});

    RhMonthsDIN20Max = maxunith(RhMonthsDIN20, minValues, maxValuesT, @(e) idealTemperatureDIN(e,20));
    hMonthsDIN20Max = maxunith(hMonthsDIN20, maxValues, minValuesT, @(e) idealTemperatureDIN(e,20));
    hMonthsDIN22Max = maxunith(hMonthsDIN22, maxValues, minValuesT, @(e) idealTemperatureDIN(e,22));
    hMonthsDIN24Max = maxunith(hMonthsDIN24, maxValues, minValuesT, @(e) idealTemperatureDIN(e,24));

    varNames = {'month', 'usage time [h]', ...
        'h < lower acceptance [h]', 'h < lower acceptance [%]', 'rev. lower acceptance', 'max. rev. lower acceptance', 'rev. lower acceptance squared', ...
        'h >= lower acceptance [h]', 'h >= lower acceptance [%]', 'lower acceptance', 'max. lower acceptance', 'lower acceptance squared', ...
        'h >= average comfort [h]', 'h >= average comfort [%]', 'average comfort', 'max. average comfort', 'average comfort squared', ...
        'h >= upper acceptance [h]', 'h >= upper acceptance [%]', 'upper acceptance', 'max. upper acceptance', 'upper acceptance squared'};
    T = table(m.', usageTimeMonths, ...
        RhMonthsDIN20, toPerc(usageTimeMonths, RhMonthsDIN20), RunithMonthsDIN20, RhMonthsDIN20Max, Runith2MonthsDIN20, ...
        hMonthsDIN20, toPerc(usageTimeMonths, hMonthsDIN20), unithMonthsDIN20, hMonthsDIN20Max, unith2MonthsDIN20, ...
        hMonthsDIN22, toPerc(usageTimeMonths, hMonthsDIN22), unithMonthsDIN22, hMonthsDIN22Max, unith2MonthsDIN22, ...
        hMonthsDIN24, toPerc(usageTimeMonths, hMonthsDIN24), unithMonthsDIN24, hMonthsDIN24Max, unith2MonthsDIN24, ...
        'VariableNames', varNames);
    writetable(T, exportFilename, 'Sheet', 'DIN EN 15251 2012-12');
    
    RhMonthsASHRAE143Max = maxunith(RhMonthsASHRAE143, minValues, maxValuesT, @(e) idealTemperatureASHRAE(e,14.3));
    hMonthsASHRAE143Max = maxunith(hMonthsASHRAE143, maxValues, minValuesT, @(e) idealTemperatureASHRAE(e,14.3));
    hMonthsASHRAE178Max = maxunith(hMonthsASHRAE178, maxValues, minValuesT, @(e) idealTemperatureASHRAE(e,17.8));
    hMonthsASHRAE213Max = maxunith(hMonthsASHRAE213, maxValues, minValuesT, @(e) idealTemperatureASHRAE(e,21.3));
    
    varNames = {'month', 'usage time [h]', ...
        'h < lower 80% acceptance [h]', 'h < lower 80% acceptance [%]', ...
        'rev. lower 80% acceptance', 'max. rev. lower 80% acceptance', 'rev. lower 80% acceptance squared', ...
        'h >= lower 80% acceptance [h]', 'h >= lower 80% acceptance [%]', ...
        'lower 80% acceptance', 'max. lower 80% acceptance', 'lower 80% acceptance squared', ...
        'h >= average comfort [h]', 'h >= average comfort [%]', ...
        'average comfort', 'max. average comfort', 'average comfort squared', ...
        'h >= upper 80% acceptance [h]', 'h >= upper 80% acceptance [%]', ...
        'upper 80% acceptance', 'max. upper 80% acceptance', 'upper 80% acceptance squared'};
    T = table(m.', usageTimeMonths, ...
        RhMonthsASHRAE143, toPerc(usageTimeMonths, RhMonthsASHRAE143), ...
        RunithMonthsASHRAE143, RhMonthsASHRAE143Max, Runith2MonthsASHRAE143, ...
        hMonthsASHRAE143, toPerc(usageTimeMonths, hMonthsASHRAE143), ...
        unithMonthsASHRAE143, hMonthsASHRAE143Max, unith2MonthsASHRAE143, ...
        hMonthsASHRAE178, toPerc(usageTimeMonths, hMonthsASHRAE178), ...
        unithMonthsASHRAE178, hMonthsASHRAE178Max, unith2MonthsASHRAE178, ...
        hMonthsASHRAE213, toPerc(usageTimeMonths, hMonthsASHRAE213), ...
        unithMonthsASHRAE213, hMonthsASHRAE213Max, unith2MonthsASHRAE213, ...
        'VariableNames', varNames);
    writetable(T, exportFilename, 'Sheet', 'ASHRAE 55-2020');
    
else
    unithMonthsMax = maxunith(hMonths, maxValues, zeros(12,1) + threshold, @(e) e);
    varNames = {'month', ...
        'usage time [h]', 'h < threshold [h]', 'h < threshold [%]', ...
        'h >= threshold [h]', 'h >= threshold [%]', ...
        'delta CO2 hours', 'max. delta CO2 hours', 'delta2 CO2 hours'};
    T = table(m.', ...
        usageTimeMonths, usageTimeMonths - hMonths, toPerc(usageTimeMonths, usageTimeMonths - hMonths), ...
        hMonths, toPerc(usageTimeMonths, hMonths), ...
        unithMonths, unithMonthsMax, unith2Months, ...
        'VariableNames', varNames);
    writetable(T, exportFilename, 'Sheet', 'CO2 hours');
end
end

function convertedDate = convertDate(input, format)
convertedDate = datenum(datetime(input, 'InputFormat', format, 'TimeZone', 'Europe/Zurich'));
end

function maxVal = maxunith(time, val, env, func)
maxVal = time .* abs(val - arrayfun(@(e) func(e), env));
end

function percentage = toPerc(total, part)
percentage = (100 * part) ./ total;
end

function filePath = selectFile(extension, title)
[file, path] = uigetfile(extension, title);
if isequal(file, 0)
    disp('Canceled');
else
    filePath = fullfile(path, file);
    disp(['User selected ', filePath]);
end
end

function selection = selectOption(options, description)
[indx, tf] = listdlg('PromptString', {description}, 'SelectionMode', 'single', 'ListString', options);
if tf == 0
    disp('Canceled');
    return;
end
selection = options(indx);
end

function xlabels = createXLabels(edges)
edges_str = string(edges);
xlabels=strcat("<", edges_str(2));
if (length(edges) > 3)
    for i = 2:(length(edges)-2)
        xlabels = [xlabels strcat(edges_str(i), "-", edges_str(i+1))];
    end
end
xlabels = [xlabels strcat(">", edges_str(length(edges) - 1))];
xlabels = reordercats(categorical(xlabels), xlabels);
end

function plotBarPlot(hist, xlabels, unit, type)
bar(xlabels, hist);
xlabel(unit);
ylabel(strcat(type, ' hours'));
text(1:length(hist),hist,strcat(num2str(hist'), " - ", num2str((hist.'/sum(hist))*100), "%"),'vert','bottom','horiz','center');
end

function filteredData = filterData(data, simulationStart, startDate, endDate, MoFr, Saturday, Sunday, holidays, DateInputFormat)
if(~strcmp(simulationStart, startDate))
    data = data(startDateInternal <= data.dates, :);
end

if(~strcmp(simulationStart, endDate))
    data = data(endDateInternal > data.dates, :);
end

sel = false(height(data),1);
if (size(MoFr) ~= 0)
    for iter = 1:height(MoFr)
        sel = sel | (1 < data.weekdayNum & data.weekdayNum < 7 & MoFr(iter,1) <= data.Hour & data.Hour < MoFr(iter,2));
    end
end

if (size(Saturday) ~= 0)
    for iter = 1:height(Saturday)
        sel = sel | (data.weekdayNum == 7 & Saturday(iter,1) <= data.Hour & data.Hour < Saturday(iter,2));
    end
end

if (size(Sunday) ~= 0)
    for iter = 1:height(Sunday)
        sel = sel | (data.weekdayNum == 1 & Sunday(iter,1) <= data.Hour & data.Hour < Sunday(iter,2));
    end
end

if (size(holidays) ~= 0)
    for iter = 1:height(holidays)
        beginDate2 = convertDate(holidays(iter, 1),DateInputFormat);
        endDate2 = convertDate(holidays(iter, 2),DateInputFormat);
        sel = sel & (data.dates < beginDate2 | data.dates > endDate2);
    end
end

filteredData = data(sel, :);
end

function [hMonths, unithMonths, unith2Months] = calculateUnitHours(data, sampleInterval, unit)
hMonths = zeros(12,1);
unithMonths = zeros(12,1);
unith2Months = zeros(12,1);
for m = 1:12
    sel = data.values >= data.threshold & data.Month == m;
    sel_size = height(data.values(sel));
    if(sel_size ~= 0)
        hMonths(m) = sel_size / (60 / sampleInterval);
        unithMonths(m) = sum(data.values(sel) - data.threshold(sel)) / (60 / sampleInterval);
        unith2Months(m) = sum((data.values(sel) - data.threshold(sel)).^2) / (60 / sampleInterval);
    end
    disp(int2str(m));
    disp(strcat("#hours: ", int2str(hMonths(m))));
    disp(strcat("Delta ", unit, " hours: ",int2str(unithMonths(m))));
    disp(strcat("Delta² ", unit, " hours: ", int2str(unith2Months(m))));
end

disp(" ");
disp(strcat("total #hours: ", int2str(sum(hMonths))));
disp(strcat("total delta ", unit, " hours: ",int2str(sum(unithMonths))));
disp(strcat("total delta² ", unit, " hours: ", int2str(sum(unith2Months))));
end

function [hMonths, unithMonths, unith2Months] = calculateUnitHoursDynThreshold(data, sampleInterval, unit, func)
hMonths = zeros(12,1);
unithMonths = zeros(12,1);
unith2Months = zeros(12,1);

dynThreshold = arrayfun(@(e) func(e),data.threshold);
for m = 1:12
    sel = data.values >= dynThreshold & data.Month == m;
    sel_size = height(data.values(sel));
    if(sel_size ~= 0)
        hMonths(m) = sel_size / (60 / sampleInterval);
        unithMonths(m) = sum(data.values(sel) - dynThreshold(sel)) / (60 / sampleInterval);
        unith2Months(m) = sum((data.values(sel) - dynThreshold(sel)).^2) / (60 / sampleInterval);
    end
end
end

function [hMonths, unithMonths, unith2Months] = calculateRevUnitHoursDynThreshold(data, sampleInterval, unit, func)
hMonths = zeros(12,1);
unithMonths = zeros(12,1);
unith2Months = zeros(12,1);

dynThreshold = arrayfun(@(e) func(e),data.threshold);
for m = 1:12
    sel = data.values < dynThreshold & data.Month == m;
    sel_size = height(data.values(sel));
    if(sel_size ~= 0)
        hMonths(m) = sel_size / (60 / sampleInterval);
        unithMonths(m) = sum(dynThreshold(sel) - data.values(sel)) / (60 / sampleInterval);
        unith2Months(m) = sum((dynThreshold(sel) - data.values(sel)).^2) / (60 / sampleInterval);
    end
end
end

% DIN EN 15251:2012-12
% lower acceptance base: 20 [°C]
% comfort base: 22 [°C]
% upper acceptance base: 24 [°C]
function temperature = idealTemperatureDIN(outsideAirTemperature, base)
if outsideAirTemperature < 16
    temperature = base;
elseif outsideAirTemperature < 32
    temperature = base + (outsideAirTemperature - 16) * 0.25;
else
    temperature = base + 4;
end
end

% ASHRAE Standard 55-2020
% https://www.ashrae.org/technical-resources/standards-and-guidelines/read-only-versions-of-ashrae-standards
% lower 80% acceptance base: 14.3 [°C]
% avg. comfort base: 17.8 [°C]
% upper 80% acceptance base: 21.3 [°C]
function temperature = idealTemperatureASHRAE(outsideAirTemperature, base)
temperature = 0.31 * outsideAirTemperature + base;
end
