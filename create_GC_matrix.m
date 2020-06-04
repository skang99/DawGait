function [] = create_GC_matrix(new_RGT_x,new_RGT_y,new_RGT_z,new_RLE_x,new_RLE_y,new_RLE_z,new_R_Tricep_x,new_R_Tricep_y,new_R_Tricep_z)

if exist('RGT_x_GC','var') == 0
     
        RGT_x_GC = new_RGT_x;        
        RGT_y_GC = new_RGT_y;        
        RGT_z_GC = new_RGT_z;
        
        RLE_x_GC = new_RLE_x;        
        RLE_y_GC = new_RLE_y;
        RLE_z_GC = new_RLE_z;
        
        R_Tricep_x_GC = new_R_Tricep_x;        
        R_Tricep_y_GC = new_R_Tricep_y;
        R_Tricep_z_GC = new_R_Tricep_z;
        
        save GaitCycles.mat RGT_x_GC RGT_y_GC RGT_z_GC RLE_x_GC RLE_y_GC RLE_z_GC R_Tricep_x_GC R_Tricep_y_GC R_Tricep_z_GC -v7.3;
       
else
    
	RGT_x_2 = AA.RGT_x_GC;        
	RGT_y_2 = AA.RGT_y_GC;        
	RGT_z_2 = AA.RGT_z_GC;
        
	RLE_x_2 = AA.RLE_x_GC;       
    RLE_y_2 = AA.RLE_y_GC;       
    RLE_z_2 = AA.RLE_z_GC;
        
    R_Tricep_x_2 = AA.R_Tricep_x_GC;        
    R_Tricep_y_2 = AA.R_Tricep_y_GC;       
    R_Tricep_z_2 = AA.R_Tricep_z_GC;
        
    t1 = time(1:length(new_RGT_x));
        
    inc = t1(end)/(length(RGT_x_2) - 1); %Finds the increment we need
                                             %for the new time array t2
        
    t2 = 0:inc:t1(end); %Creates an array of the new time that will be
                            %used to interpolate between the two columns of data
        
        
    new_RGT_x_2 = interp1(t1,new_RGT_x,t2);       
    new_RGT_y_2 = interp1(t1,new_RGT_y,t2);       
    new_RGT_z_2 = interp1(t1,new_RGT_z,t2);
        
    new_RLE_x_2 = interp1(t1,new_RLE_x,t2);        
    new_RLE_y_2 = interp1(t1,new_RLE_y,t2);       
    new_RLE_z_2 = interp1(t1,new_RLE_z,t2);
        
    new_R_Tricep_x_2 = interp1(t1,new_R_Tricep_x,t2);       
    new_R_Tricep_y_2 = interp1(t1,new_R_Tricep_y,t2);        
    new_R_Tricep_z_2 = interp1(t1,new_R_Tricep_z,t2);
        
    %The interpolation function converts the position matrices from
    %column vectors into row vectors. In order to keep everything
    %consistent, the following code transforms the matrices back into
    %column vectors
        
    new_RGT_x = new_RGT_x_2.';       
    new_RGT_y = new_RGT_y_2.';        
    new_RGT_z = new_RGT_z_2.';
        
    new_RLE_x = new_RLE_x_2.';        
    new_RLE_y = new_RLE_y_2.';
    new_RLE_z = new_RLE_z_2.';
        
    new_R_Tricep_x = new_R_Tricep_x_2.';        
    new_R_Tricep_y = new_R_Tricep_y_2.';        
    new_R_Tricep_z = new_R_Tricep_z_2.';
           
    RGT_x_GC(:,size(RGT_x_2,2) + 1) = new_RGT_x;
    RGT_y_GC(:,size(RGT_x_2,2) + 1) = new_RGT_y;       
    RGT_z_GC(:,size(RGT_x_2,2) + 1) = new_RGT_z;
        
    RLE_x_GC(:,size(RGT_x_2,2) + 1) = new_RLE_x;       
    RLE_y_GC(:,size(RGT_x_2,2) + 1) = new_RLE_y;        
    RLE_z_GC(:,size(RGT_x_2,2) + 1) = new_RLE_z;
        
    R_Tricep_x_GC(:,size(R_Tricep_x_2,2) + 1) = new_R_Tricep_x;        
    R_Tricep_y_GC(:,size(R_Tricep_y_2,2) + 1) = new_R_Tricep_y;       
    R_Tricep_z_GC(:,size(R_Tricep_z_2,2) + 1) = new_R_Tricep_z;

    save GaitCycles.mat RGT_x_GC RGT_y_GC RGT_z_GC RLE_x_GC RLE_y_GC RLE_z_GC R_Tricep_x_GC R_Tricep_y_GC R_Tricep_z_GC  -v7.3;    
    
end   
    
end

