function [static_RGT_RLE,static_RLO_RLS,static_RSC1_RAC,static_RDS_RAC,static_ACB_R5M,T1,VTR1,SC1,DLMC5,ACB,RDS,RTR,RME,RCR,R5M,R2M,RMS,RGT,RLE,Centroid,RLS,RLO,VTR2,VTR3,RCDS] = create_static_data(static_trial)

% Extracts positional data from the xls file static_trial
% Assumes xls files passed in are all similiarly formatted, namely,
% that the row data begins in is row 6
% static_trial should not containing leading blank/0 data

ROW_START = 6;

try 
[static_data,static_coords] = xlsread(static_trial); 
catch exception
    disp("Error in opening the specified static Excel file. File may not exist or is not within the current directory.");
    return
end

%New system
static_pos_data = static_data(ROW_START:length(static_data),:); 

RGT = double(subs(extract_data(static_pos_data,static_coords,'GRTB'),NaN,0));
RLE = double(subs(extract_data(static_pos_data,static_coords,'LEPI'),NaN,0));
R5M = double(subs(extract_data(static_pos_data,static_coords,'MCP5'),NaN,0));
R2M = double(subs(extract_data(static_pos_data,static_coords,'MCP2'),NaN,0));
RLO = double(subs(extract_data(static_pos_data,static_coords,'OLEC'),NaN,0));
RLS = double(subs(extract_data(static_pos_data,static_coords,'LSTY'),NaN,0));
T1 = double(subs(extract_data(static_pos_data,static_coords,'T1'),NaN,0));
RDS = double(subs(extract_data(static_pos_data,static_coords,'DRSC'),NaN,0));
RME = double(subs(extract_data(static_pos_data,static_coords,'MEPI'),NaN,0));
RTR = double(subs(extract_data(static_pos_data,static_coords,'TRIC'),NaN,0));
RCR = double(subs(extract_data(static_pos_data,static_coords,'MDCR'),NaN,0));
RMS = double(subs(extract_data(static_pos_data,static_coords,'MSTY'),NaN,0));

RAC = double(subs(extract_data(static_pos_data,static_coords,'ACRM'),NaN,0));
RSC1 = double(subs(extract_data(static_pos_data,static_coords,'SC1'),NaN,0));
RSC2 = double(subs(extract_data(static_pos_data,static_coords,'SC2'),NaN,0));

VTR1 = double(subs(extract_data(static_pos_data,static_coords,'VTR1'),NaN,0));
VTR2 = double(subs(extract_data(static_pos_data,static_coords,'VT2'),NaN,0));
VTR3 = double(subs(extract_data(static_pos_data,static_coords,'TV3'),NaN,0));
RCDS = double(subs(extract_data(static_pos_data,static_coords,'CDSC'),NaN,0));

SC1 = double(subs(extract_data(static_pos_data,static_coords,'SC1'),NaN,0));
DLMC5 = double(subs(extract_data(static_pos_data,static_coords,'DLMC5'),NaN,0));
ACB = double(subs(extract_data(static_pos_data,static_coords,'ACCB'),NaN,0));

Centroid = (RAC + RSC1 + RSC2) / 3;

lower_lim =  1; 
upper_lim =  100;

% Assumes the first 100 frames are when the static readings vary the least
% and that data is not blank

T1 = T1(lower_lim: upper_lim,:);
VTR1 = VTR1(lower_lim: upper_lim,:);
SC1 = SC1(lower_lim: upper_lim,:);
DLMC5 = DLMC5(lower_lim: upper_lim,:);
ACB = ACB(lower_lim: upper_lim,:);
RDS = RDS(lower_lim: upper_lim,:);
RAC = RAC(lower_lim: upper_lim,:);
RTR = RTR(lower_lim: upper_lim,:);
RME = RME(lower_lim: upper_lim,:);
RCR = RCR(lower_lim: upper_lim,:);
R5M = R5M(lower_lim: upper_lim,:);
R2M = R2M(lower_lim: upper_lim,:);
RMS = RMS(lower_lim: upper_lim,:);
RGT = RGT(lower_lim: upper_lim,:);
RLE = RLE(lower_lim: upper_lim,:);
Centroid = Centroid(lower_lim: upper_lim,:);
RLS = RLS(lower_lim: upper_lim,:);
RLO = RLO(lower_lim: upper_lim,:);
VTR2 = VTR2(lower_lim: upper_lim,:);
VTR3 = VTR3(lower_lim: upper_lim,:); 
RCDS = RCDS(lower_lim: upper_lim,:);
RSC1 = RSC1(lower_lim: upper_lim,:);


ab_x = RGT(:,1) - RLE(:,1); 
ab_z = RGT(:,3) - RLE(:,3); 
ab = sqrt(ab_x.^2 + ab_z.^2);

%Averages the previous distances to create a scalar. 
%This scalar value will be used as a baseline to
%determine whether or not dynamic gait cycles
%have too much skin movement in the error_check function.

static_RGT_RLE = mean(ab); 
                      
de_x = RLO(:,1) - RLS(:,1);
de_z = RLO(:,3) - RLS(:,3);
de = sqrt(de_x.^2 + de_z.^2); 
static_RLO_RLS = mean(de); 

gj_x = RSC1(:,1) - RAC(:,1); 
gj_z = RSC1(:,3) - RAC(:,3); 
gj = sqrt(gj_x.^2 + gj_z.^2); 
static_RSC1_RAC = mean(gj);                       
                           
hb_x = RDS(:,1) - RAC(:,1); 
hb_z = RDS(:,3) - RAC(:,3); 
hb = sqrt(hb_x.^2 + hb_z.^2);
static_RDS_RAC = mean(hb);

am_x = ACB(:,1) - R5M(:,1);
am_z = ACB(:,3) - R5M(:,3);
am = sqrt(am_x.^2 + am_z.^2);
static_ACB_R5M = mean(am);


                                                                              
end







