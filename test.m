
gc_count = 3;

filename = "test angular Data.xlsx";

name = {'Dog 1'
        ''
        'GC#'};

header = {'1' '' '' '' '' '' '' '' ''};
label = {'ELB_f' 'ELB_r' 'ELB_a' 'CARP_f' 'CARP_r' 'CARP_a' 'SHLD_f' 'SHLD_r' 'SHLD_a'};

header = repmat(header,1,gc_count);
label = repmat(label,1,gc_count);
index = 10;

for i=1:gc_count - 1
    header(index) = num2cell(i + 1);
    index = index + 9;
end

writecell(name,filename,'Range','A1')
writecell(header,filename,'Range','B3');
writecell(label,filename,'Range','B4');






