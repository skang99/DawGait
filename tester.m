% 
% t1  = [1;2;3;4;5;6;7;8;9];
% t2 = [10;11;12;13;14;15;16;17;18;19];
% 
% name = {'Colt No Vest 04'
%         ''
%         'GC#'};
% header = {'1' '' '' '' '2' '' '' '' '3'
%            'BA_f' 'BA_r' 'BA_a' 'AP_f' 'AP_r' 'AP_a' 'BS_f' 'BS_r' 'BS_a'};
% 
% filename = 'dogtest.xlsx';
% T = table(jc_angles.angles(1).ba_f,jc_angles.angles(1).bs_f,jc_angles.angles(1).ap_f,jc_angles.angles(1).ba_r,jc_angles.angles(1).bs_r,jc_angles.angles(1).ap_r,jc_angles.angles(1).ba_a,jc_angles.angles(1).bs_a,jc_angles.angles(1).ap_a)
% 
% 
% %Range arguments possibly accept variables?
% %Could be able to append by incrementing this value
% writecell(name,filename,'Range','A1')
% writecell(header,filename,'Range','B3')
% 
% writetable(T,filename,'Sheet',1,'Range','B5','WriteVariableNames',false)

% x = [1:10]
% 
% save('Produced Data/tester.mat','x')

% 

x = [1:30];
y = [1:50];

figure(1)
plot(x,'Color','#f55151')

figure(2)
plot(y,'Color','#5192f5')




