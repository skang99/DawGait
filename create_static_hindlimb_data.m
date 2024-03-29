function [static_RGT_RLEP,static_RFH_RLMA,static_RIWG_RISC,static_CRS_RISC,static_RCAL_RMP5,RGT,RLEP,RMP5,MP2,RFH,RLMA,CRS,RIWG,RMEP,RQUA,RGAS,RISC,RMMA,RCAL,DLMC5,CGT,PTC] = create_static_hindlimb_data(static_trial)

% Extracts positional data from the xls file static_trial
% Assumes xls files passed in are all similiarly formatted, namely,
% that the row data begins in is row 6
% static_trial should not containing leading blank/0 data

ROW_START = 6;

try 
[static_data,static_coords] = xlsread(static_trial); %import excel data numeric and text
catch exception
    disp("Error in opening the specified static Excel file. File may not exist or is not within the current directory.");
    return;
end

%New system
static_pos_data = static_data(ROW_START:length(static_data),:); %Extract only the numeric data a.k.a position readings


RGT = double(subs(extract_data(static_pos_data,static_coords,'GT'),NaN,0));
RLEP = double(subs(extract_data(static_pos_data,static_coords,'LEP'),NaN,0));
RMP5 = double(subs(extract_data(static_pos_data,static_coords,'MP5'),NaN,0));
RFH = double(subs(extract_data(static_pos_data,static_coords,'FH'),NaN,0));
RLMA = double(subs(extract_data(static_pos_data,static_coords,'LMA'),NaN,0));
CRS = double(subs(extract_data(static_pos_data,static_coords,'RS'),NaN,0));
RIWG = double(subs(extract_data(static_pos_data,static_coords,'IWG'),NaN,0));
RMEP = double(subs(extract_data(static_pos_data,static_coords,'MEP'),NaN,0));
RQUA = double(subs(extract_data(static_pos_data,static_coords,'QUA'),NaN,0));
RGAS = double(subs(extract_data(static_pos_data,static_coords,'GAS'),NaN,0));
RISC = double(subs(extract_data(static_pos_data,static_coords,'ISC'),NaN,0));
RMMA = double(subs(extract_data(static_pos_data,static_coords,'MMA'),NaN,0));
RCAL = double(subs(extract_data(static_pos_data,static_coords,'CAL'),NaN,0));
CGT =double(subs(extract_data(static_pos_data,static_coords,'CGT'),NaN,0));
PTC = double(subs(extract_data(static_pos_data,static_coords,'PTC'),NaN,0));
DLMC5 = double(subs(extract_data(static_pos_data,static_coords,'DLMC5'),NaN,0));
MP2 = double(subs(extract_data(static_pos_data,static_coords,'MP2'),NaN,0));

% Assumes the first 100 frames are when the static readings vary the least

lower_lim =  1; 
upper_lim =  100;

RGT = RGT(lower_lim: upper_lim,:); 
RLEP = RLEP(lower_lim: upper_lim,:); 
RMP5 = RMP5(lower_lim: upper_lim,:);
RFH = RFH(lower_lim: upper_lim,:); 
RLMA = RLMA(lower_lim: upper_lim,:); 
RIWG = RIWG(lower_lim: upper_lim,:);
RCAL = RCAL(lower_lim:upper_lim,:);
RISC = RISC(lower_lim:upper_lim,:);
CRS = CRS(lower_lim:upper_lim,:);

ab_x = RGT(:,1) - RLEP(:,1); 
ab_z = RGT(:,3) - RLEP(:,3); 
ab = sqrt(ab_x.^2 + ab_z.^2);

%Averages the previous distances to create a scalar. 
%This scalar value will be used as a baseline to
%determine whether or not dynamic gait cycles
%have too much skin movement in the error_check function.
static_RGT_RLEP = mean(ab);
                      
de_x = RFH(:,1) - RLMA(:,1);
de_z = RFH(:,3) - RLMA(:,3);
de = sqrt(de_x.^2 + de_z.^2); 
static_RFH_RLMA = mean(de); 

gj_x = RIWG(:,1) - RISC(:,1); 
gj_z = RIWG(:,3) - RISC(:,3); 
gj = sqrt(gj_x.^2 + gj_z.^2); 
static_RIWG_RISC = mean(gj);      
              
hb_x = CRS(:,1) - RISC(:,1); 
hb_z = CRS(:,3) - RISC(:,3); 
hb = sqrt(hb_x.^2 + hb_z.^2);
static_CRS_RISC = mean(hb);

am_x = RCAL(:,1) - RMP5(:,1);
am_z = RCAL(:,3) - RMP5(:,3);
am = sqrt(am_x.^2 + am_z.^2);
static_RCAL_RMP5 = mean(am);

                                                                                
end







