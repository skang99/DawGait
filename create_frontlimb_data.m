function [R5M,R2M,RGT,RLE,RLO,RLS,T1,RDS,Centroid,time,RME,RMS,RTR,RCR,ACB,RAC,DLMC5,VTR1,VTR2,VTR3,RSC1,RCDS] = create_frontlimb_data(filename)
% Extracts positional data from the xls file filename
% Assumes xls files passed in are all similiarly formatted, namely,
% that the row data begins in is row 6

ROW_START = 6;

try 
    [kinematic_data,markers] = xlsread(filename);
catch exception
    disp("Error in opening the specified Excel file. File may not exist or is not within the current directory.");
    disp(exception)
    throw(exception)
end

%Creates vectors for time and position data 
position_data = kinematic_data(ROW_START:length(kinematic_data),:);
time = position_data(:,1);

time = time / 200;

RGT = double(subs(extract_data(position_data,markers,'GRTB'),NaN,0));
R2M = double(subs(extract_data(position_data,markers,'MCP2'),NaN,0));
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
DLMC5 = double(subs(extract_data(position_data,markers,'DLMC5'),NaN,0));

VTR1 = double(subs(extract_data(position_data,markers,'VTR1'),NaN,0));
VTR2 = double(subs(extract_data(position_data,markers,'VT2'),NaN,0));
VTR3 = double(subs(extract_data(position_data,markers,'TV3'),NaN,0));
RCDS = double(subs(extract_data(position_data,markers,'CDSC'),NaN,0));


Centroid = (RAC + RSC1 + RSC2) * (1/3);



end

