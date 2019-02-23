% read in .rpo files from oasn and generate the replicas from
% each file in the directory

array_size = 15;
prefix = '/Users/Rui/Desktop/Ice_source/1720_int0.00008/';
directory = dir([prefix '*.txt']);

oasn_rpo = zeros(array_size,length(directory));
rpo_cov = zeros(array_size,array_size,length(directory));

for k = 1:length(directory)
    
    filename = [prefix directory(k).name];
    fileID = fopen(filename);
    C = textscan(fileID,'%s%s%s%s%s%s%s%s%s','HeaderLines',array_size+14);

    for i = 1:array_size
        
        re(i) = str2double(C{1,7}(i));
        im(i) = str2double(C{1,8}{i,1}(1:end-1));
    end
    
    oasn_rpo(:,k) = complex(re,im);
    rpo_cov(:,:,k) = oasn_rpo(:,k)*oasn_rpo(:,k)';

end

