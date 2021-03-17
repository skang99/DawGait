function [rMarker] = reconstruct_marker(m1,m2,m3,sm1,sm2,sm3,sm4)
% sm4 should be the static data for the corresponding rMarker
% all other dynamic (m) and static (sm) markers should correspond


[m1_x m1_y m1_z] = extract_XYZ(m1);
[m2_x m2_y m2_z] = extract_XYZ(m2);
[m3_x m3_y m3_z] = extract_XYZ(m3);
    
[sm1_x sm1_y sm1_z] = extract_XYZ(sm1);
[sm2_x sm2_y sm2_z] = extract_XYZ(sm2);
[sm3_x sm3_y sm3_z] = extract_XYZ(sm3);
[sm4_x sm4_y sm4_z] = extract_XYZ(sm4);

sm1_x = mean(sm1_x); sm1_y = mean(sm1_y); sm1_z = mean(sm1_z);
sm2_x = mean(sm2_x); sm2_y = mean(sm2_y); sm2_z = mean(sm2_z);
sm3_x = mean(sm3_x); sm3_y = mean(sm3_y); sm3_z = mean(sm3_z);
sm4_x = mean(sm4_x); sm4_y = mean(sm4_y); sm4_z = mean(sm4_z);

%temp reconstructed marker
trm1_x = m1_x - sm1_x + sm4_x;
trm1_y = m1_y - sm1_y + sm4_y;
trm1_z = m1_z - sm1_z + sm4_z;

rm1 = [trm1_x trm1_y trm1_z];

trm2_x = m2_x - sm2_x + sm4_x;
trm2_y = m2_y - sm2_y + sm4_y;
trm2_z = m2_z - sm2_z + sm4_z;

rm2 = [trm2_x trm2_y trm2_z];

trm3_x = m3_x - sm3_x + sm4_x;
trm3_y = m3_y - sm3_y + sm4_y;
trm3_z = m3_z - sm3_z + sm4_z;

rm3 = [trm3_x trm3_y trm3_z];

rMarker = [((trm1_x + trm2_x + trm3_x) / 3) ((trm2_x + trm2_x + trm3_x) / 3) ((trm3_x + trm2_x + trm3_x) / 3)];

end

