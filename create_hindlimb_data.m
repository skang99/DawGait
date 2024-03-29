function [RMP5,RMP2,RGT,RLEP,RFH,RLMA,CRS,RIWG,time,RMEP,RMMA,RQUA,RGAS,RCAL,RISC,CGT,PTC,side] = create_hindlimb_cycles(filename)
% Extracts positional data from the xls file filename
% Assumes xls files passed in are all similiarly formatted, namely,
% that the row data begins in is row 6

ROW_START = 6;

try 
    [kinematic_data,markers] = xlsread(filename);
catch exception
    disp("Error in opening the specified Excel file. File may not exist or is not within the current directory.");
    disp(exception)
    throw exception
end

% Addresses 3rd column, 3rd row of spreadsheet (Cell C3)
if(cell2mat(strfind(markers(3,3),"LIWG")))
    side = -1; %left
else
    side = 1; %right
end
    

%Creates vectors for time and position data 
position_data = kinematic_data(ROW_START:length(kinematic_data),:);
time = position_data(:,1);

time = time / 200;

RMP2 = double(subs(extract_data(position_data,markers,'MP2'),NaN,0));
RGT = double(subs(extract_data(position_data,markers,'GT'),NaN,0));
RLEP = double(subs(extract_data(position_data,markers,'LEP'),NaN,0));
RMP5 = double(subs(extract_data(position_data,markers,'MP5'),NaN,0));
RFH = double(subs(extract_data(position_data,markers,'FH'),NaN,0));
RLMA = double(subs(extract_data(position_data,markers,'LMA'),NaN,0));
CRS = double(subs(extract_data(position_data,markers,'CRS'),NaN,0));
RIWG = double(subs(extract_data(position_data,markers,'IWG'),NaN,0));
RMEP = double(subs(extract_data(position_data,markers,'MEP'),NaN,0));
RQUA = double(subs(extract_data(position_data,markers,'QUA'),NaN,0));
RGAS = double(subs(extract_data(position_data,markers,'GAS'),NaN,0));

RISC = double(subs(extract_data(position_data,markers,'ISC'),NaN,0));

RMMA = double(subs(extract_data(position_data,markers,'MMA'),NaN,0));
RCAL = double(subs(extract_data(position_data,markers,'CAL'),NaN,0));
CGT = double(subs(extract_data(position_data,markers,'CGT'),NaN,0));
PTC = double(subs(extract_data(position_data,markers,'PTC'),NaN,0));
 

end

