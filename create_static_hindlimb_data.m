function [static_RGT_RLEP,static_RFH_RLMA,static_RIWG_RISC,static_CRS_RISC,static_RCAL_RMP5,RGT,RLEP,RMP5,RFH,RLMA,CRS,RIWG,RMEP,RQUA,RGAS,RISC,RMMA,RCAL,DLMC5] = create_static_hindlimb_data(static_trial,new_new_time,time)

try 
[static_data,static_coords] = xlsread(static_trial); %import excel data numeric and text
catch exception
    disp("Error in opening the specified static Excel file. File may not exist or is not within the current directory.");
    return;
end

%New system
static_pos_data = static_data(6:length(static_data),:); %Extract only the numeric data a.k.a position readings


RGT = double(subs(extract_data(static_pos_data,static_coords,'RGT'),NaN,0));
RLEP = double(subs(extract_data(static_pos_data,static_coords,'RLEP'),NaN,0));
RMP5 = double(subs(extract_data(static_pos_data,static_coords,'RMP5'),NaN,0));
RFH = double(subs(extract_data(static_pos_data,static_coords,'RFH'),NaN,0));
RLMA = double(subs(extract_data(static_pos_data,static_coords,'RLMA'),NaN,0));
CRS = double(subs(extract_data(static_pos_data,static_coords,'CRS'),NaN,0));
RIWG = double(subs(extract_data(static_pos_data,static_coords,'RIWG'),NaN,0));
RMEP = double(subs(extract_data(static_pos_data,static_coords,'RMEP'),NaN,0));
RQUA = double(subs(extract_data(static_pos_data,static_coords,'RQUA'),NaN,0));
RGAS = double(subs(extract_data(static_pos_data,static_coords,'RGAS'),NaN,0));

RISC = double(subs(extract_data(static_pos_data,static_coords,'RISC'),NaN,0));


RMMA = double(subs(extract_data(static_pos_data,static_coords,'RMMA'),NaN,0));
RCAL = double(subs(extract_data(static_pos_data,static_coords,'RCAL'),NaN,0));

DLMC5 = double(subs(extract_data(static_pos_data,static_coords,'DLMC5'),NaN,0));


R_5th_M_z = RMP5(:,3);

[t_naught, t_one] = find_still_static(R_5th_M_z);


% lower_lim =  %finds the frame where the lower limit time value is located
% upper_lim = find(time == t_one) %finds the frame where the upper limit time value is located

lower_lim =  1; %finds the frame where the lower limit time value is located
upper_lim =  100;%finds the frame where the upper limit time value is located

RGT = RGT(lower_lim: upper_lim,:); %extracts each frame of data within the specified range
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







