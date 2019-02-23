% read in .chk files from oasn and generate the covariance matrices from
% each file in the directory

array_size = 15;
prefix = '/Users/Rui/Desktop/Ice_source/chk_files/';
directory = dir([prefix '*.chk']);

oasn_cov = zeros(array_size,array_size,length(directory));

for k = 1:length(directory)

    filename = [prefix directory(k).name];
    fileID = fopen(filename);
    C = textscan(fileID,'%*s%s%s%s%*s%*s%*s','HeaderLines',array_size+9);

    for i = 1:length(C{1,1});
        num(i) = str2double(C{1,1}(i));
        re(i) = str2double(C{1,2}(i));
        im(i) = str2double(C{1,3}(i));
    end

    keep = ~isnan(num) & (floor(num) == num);
    re = re(keep);
    im = im(keep);

    data = complex(re,im);

    counter = 1;
    for i = 1:array_size
        for j = i:array_size
            oasn_cov(i,j,k) = data(counter);

            if i ~= j
                oasn_cov(j,i,k) = conj(data(counter));
            end

            counter = counter + 1;
        end
    end
end


num_int = 10;

oasn_covb1 = oasn_cov;
oasn_covb2 = zeros(array_size,array_size,length(directory));

counter = 0;
for int = 1:num_int
    count = 0;
    for int2 = 1:num_int
        oasn_covb2(:,:,int2+counter*num_int) = oasn_cov(:,:,int+count);
        count = count+num_int;
    end
    count = count+1;
    counter = counter+1;
    
end
    
    