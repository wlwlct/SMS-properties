%use to test whether the dtime start at same place
clearvars
cd('E:\F8T2400nmCH\apd full')
allnames=struct2cell(dir('F8T2*'));
[~,len]=size(allnames);

total_SecDtime02=[];
total_SecDtime03=[];
total_SecDtime04=[];
total_SecDtime05=[];
for len_i=1:1:len
    date=regexp(allnames{1,len_i},'02\d*2020','match');
    
    switch date{1}
    
        case '02022020'
            total_SecDtime02=[total_SecDtime02;importdata(char(allnames(1,len_i)))];
        case '02032020'
            total_SecDtime03=[total_SecDtime03;importdata(char(allnames(1,len_i)))];
        case '02042020'
            total_SecDtime04=[total_SecDtime04;importdata(char(allnames(1,len_i)))];
        case '02052020'
            total_SecDtime05=[total_SecDtime05;importdata(char(allnames(1,len_i)))];
        otherwise
            disp('Not match with anydate')
    end
end

test02=cell2mat(total_SecDtime02(:,2));
test03=cell2mat(total_SecDtime03(:,2));
test04=cell2mat(total_SecDtime04(:,2));
test05=cell2mat(total_SecDtime05(:,2));

%Check dtime during each day
%figure;mesh(test02(:,100:500));view([0 0 1])

%Check dtime difference during those days
test02_sum=sum(test02,1);
test03_sum=sum(test03,1);
test04_sum=sum(test04,1);
test05_sum=sum(test05,1);

figure;plot(test02_sum)
hold on;plot(test03_sum)
hold on;plot(test04_sum)
hold on;plot(test05_sum)
xlim([100 500])

figure;plot(normalize(test02_sum,2,'range'))
figure;plot(normalize(test03_sum,2,'range'))
hold on;plot(normalize(test03_sum,2,'range'))
hold on;plot(normalize(test04_sum,2,'range'))
hold on;plot(normalize(test05_sum,2,'range'))
xlim([100 500])
