function [RGT_RLE_seg,RTr_RLE_seg,RLO_RLS_seg,T1_T13_seg,RCR_RLS_seg,RAC_RDS_seg,RAC_RLE_seg,RSC1_RDS_seg,RSC2_RDS_seg,RSC1_RLE_seg,RSC2_RLE_seg] = create_dynamic_segs(RGT,R_Tricep,RLO,RLS,RCR,RLE,RAC,T1,T13,RDS,RSC1,RSC2)
%Upper Limb 

RGT_RLE_seg = sqrt((RGT(:,1) - RLE(:,1)).^2 + (RGT(:,3) - RLE(:,3)).^2);

RTr_RLE_seg = sqrt((R_Tricep(:,1) - RLE(:,1)).^2 + (R_Tricep(:,3) - RLE(:,3)).^2);

T1_T13_seg = sqrt((T1(:,1) - T13(:,1)).^2 + (T1(:,3) - T13(:,3)).^2);

%Lower Limb

RLO_RLS_seg = sqrt((RLO(:,1) - RLS(:,1)).^2 + (RLO(:,3) - RLS(:,3)).^2);

RCR_RLS_seg = sqrt((RLS(:,1) - RCR(:,1)).^2 + (RLS(:,3) - RCR(:,3)).^2);

RAC_RDS_seg = sqrt((RAC(:,1) - RDS(:,1)).^2 + (RAC(:,3) - RDS(:,3)).^2);

RAC_RLE_seg = sqrt((RAC(:,1) - RLE(:,1)).^2 + (RAC(:,3) - RLE(:,3)).^2);

RSC1_RDS_seg = sqrt((RSC1(:,1) - RDS(:,1)).^2 + (RSC1(:,3) - RDS(:,3)).^2);

RSC2_RDS_seg = sqrt((RSC2(:,1) - RDS(:,1)).^2 + (RSC2(:,3) - RDS(:,3)).^2);

RSC1_RLE_seg = sqrt((RSC1(:,1) - RLE(:,1)).^2 + (RSC1(:,3) - RLE(:,3)).^2);

RSC2_RLE_seg = sqrt((RSC2(:,1) - RLE(:,1)).^2 + (RSC2(:,3) - RLE(:,3)).^2);


end

