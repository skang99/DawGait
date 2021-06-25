
function [RMP5,RGT,RLEP,RFH,RLMA,CRS,RIWG,time,RMEP,RMMA,RQUA,RGAS,RCAL,RISC,direction] = hindlimbJCSextract(filename)
%This function takes an xls file and creates graphical
%representations of the data input. 
%filename: name of the xlsx file to read data from
%RETURNS: vector of position data for each specified landmark

%Reads the xls file and extracts position data and marker designations
%Surrounded by a try-catch block to handle potential errors
try 
    [kinematic_data,markers] = xlsread(filename);
catch exception
    rethrow(exception)
end


%Creates vectors for time and position data 
position_data = kinematic_data(6:length(kinematic_data),:);
time = position_data(:,1);

% Stores time as frames/fps 
time = time / 200;

RGT = double(subs(extract_data(position_data,markers,'GT'),NaN,0));
RLEP = double(subs(extract_data(position_data,markers,'LEP'),NaN,0));
RMP5 = double(subs(extract_data(position_data,markers,'MP5'),NaN,0));
RFH = double(subs(extract_data(position_data,markers,'FH'),NaN,0));
RLMA = double(subs(extract_data(position_data,markers,'LMA'),NaN,0));
CRS = double(subs(extract_data(position_data,markers,'RS'),NaN,0));
RIWG = double(subs(extract_data(position_data,markers,'IWG'),NaN,0));
RMEP = double(subs(extract_data(position_data,markers,'MEP'),NaN,0));
RQUA = double(subs(extract_data(position_data,markers,'QUA'),NaN,0));
RGAS = double(subs(extract_data(position_data,markers,'GAS'),NaN,0));
RISC = double(subs(extract_data(position_data,markers,'ISC'),NaN,0));
RMMA = double(subs(extract_data(position_data,markers,'MMA'),NaN,0));
RCAL = double(subs(extract_data(position_data,markers,'CAL'),NaN,0));




end

