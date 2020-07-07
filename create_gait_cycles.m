function [R5M,RGT,RLE,RLO,RLS,T1,RDS,Centroid,time,RME,RMS,RTR,RCR,ACB,RAC] = create_gait_cycles(filename)
%This function takes an xls file and creates graphical
%representations of the data input. 
%filename: name of the xlsx file to read data from
%RETURNS: vector of position data for each specified landmark

%Reads the xls file and extracts position data and marker designations
%Surrounded by a try-catch block to handle potential errors
try 
    [kinematic_data,markers] = xlsread(filename);
catch exception
    disp("Error in opening the specified Excel file. File may not exist or is not within the current directory.");
    throw exception
end

%Creates vectors for time and position data 
position_data = kinematic_data(6:length(kinematic_data),:);
time = position_data(:,1);

% position_data = kinematic_data(7:length(kinematic_data),:);
% time = position_data(:,2);



time = time / 200;

RGT = double(subs(extract_data(position_data,markers,'GRTB'),NaN,0));
RLE = double(subs(extract_data(position_data,markers,'LEPI'),NaN,0));
R5M = double(subs(extract_data(position_data,markers,'MCP5'),NaN,0));
RLO = double(subs(extract_data(position_data,markers,'OLEC'),NaN,0));
RLS = double(subs(extract_data(position_data,markers,'LSTY'),NaN,0));
T1 = double(subs(extract_data(position_data,markers,'T1'),NaN,0));
RDS = double(subs(extract_data(position_data,markers,'DRSC'),NaN,0));
RME = double(subs(extract_data(position_data,markers,'MEPI'),NaN,0));
RTR = double(subs(extract_data(position_data,markers,'TRIC'),NaN,0));
RCR = double(subs(extract_data(position_data,markers,'MDCR'),NaN,0));

RAC = double(subs(extract_data(position_data,markers,'ACRM'),NaN,0));
RSC1 = double(subs(extract_data(position_data,markers,'SC1'),NaN,0));
RSC2 = double(subs(extract_data(position_data,markers,'SC2'),NaN,0));

RMS = double(subs(extract_data(position_data,markers,'MSTY'),NaN,0));
ACB = double(subs(extract_data(position_data,markers,'ACCB'),NaN,0));

RAC = double(subs(extract_data(position_data,markers,'ACRM'),NaN,0));

% RGT = extract_data(position_data,markers,'R Greater Tubercle');
% RLE = extract_data(position_data,markers,'R Lateral Epicondyle');
% R5M = extract_data(position_data,markers,'R 5th Metacarpal');
% RLO = extract_data(position_data,markers,'R Lateral Olecranon');
% RLS = extract_data(position_data,markers,'R Lateral Styloid');
% T1 = extract_data(position_data,markers,'T1');
% RDS = extract_data(position_data,markers,'R Dorsal Scapula');
%  
% RAC = extract_data(position_data,markers,'R Acromion');
% RSC1 = extract_data(position_data,markers,'RSC1');
% RSC2 = extract_data(position_data,markers,'RSC2');
% 
% RME = extract_data(position_data,markers,'R Medial Epicondyle');
% RTR = extract_data(position_data,markers,'R Tricep');
% RMS = extract_data(position_data,markers,'R Medial Styloid');
%  


Centroid = (RAC + RSC1 + RSC2) * (1/3);



end

