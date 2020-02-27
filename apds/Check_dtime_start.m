%use to test whether the dtime start at same place
clearvars
cd('E:\F8T2N2\apd full')
allnames=struct2cell(dir('F8T2*'));
[~,len]=size(allnames);

%days=['02142019';'02152019';'02172019';'02182019';'02212019';'02222019';'02262019'];
days=['02042019'];
days_leng=length(days(:,1));

for days_i=1:days_leng
eval(['total_SecDtime' num2str(days(days_i,:)) '=[]']);
eval(['Mol_total_SecDtime' num2str(days(days_i,:)) '=zeros(len,6251);'])
end

for len_i=1:1:len
    date=regexp(allnames{1,len_i},'\d*\d*2019','match');
    for days_i=1:days_leng
        if strcmp(date{1},num2str(days(days_i,:)))
            current=importdata(char(allnames(1,len_i)));
            eval(['Mol_total_SecDtime' num2str(days(days_i,:)) '(len_i,:)=sum(cell2mat(current(:,2)),1);'])
            eval(['total_SecDtime' num2str(days(days_i,:)) '=[total_SecDtime' num2str(days(days_i,:)) '; current];'])
        end
    end
    
end

for days_i=1:days_leng
    clearvars loc
    eval(['loc=sum(Mol_total_SecDtime' num2str(days(days_i,:)) ',2)~=0;'])
    eval(['Mol_total_SecDtime' num2str(days(days_i,:)) '= Mol_total_SecDtime' num2str(days(days_i,:)) '(loc,:);'])
end

for days_i=1:days_leng
eval(['SecDtime_matrix_' num2str(days(days_i,:)) '=cell2mat(total_SecDtime' num2str(days(days_i,:)) '(:,2));']);
eval(['SecDtime_matrix_' num2str(days(days_i,:)) '_sum=sum(SecDtime_matrix_' num2str(days(days_i,:)) ',1);']);
end




%Check dtime during each sec
%figure;mesh(test02(:,100:500));view([0 0 1])
for days_i=1:days_leng
    figure;
    eval(['mesh(SecDtime_matrix_' num2str(days(days_i,:)) '(:,100:500))'])
    view([0 0 1])
end

for days_i=1:days_leng
    figure;
    eval(['mesh(normalize(SecDtime_matrix_' num2str(days(days_i,:)) '(:,100:500),2,' char(39) 'range' char(39) '))'])
    view([0 0 1])
end

%Check dtime difference based on molecule
for days_i=1:days_leng
    clearvars Mol_data
    eval(['Mol_data=Mol_total_SecDtime' days(days_i,:) ';']);
    Mol_data_leng=length(Mol_data(:,1));
    figure
    for Mol_data_i=1:Mol_data_leng
        hold on;plot(Mol_data(Mol_data_i,:));xlim([0 600])
    end
end

for days_i=1:days_leng
    clearvars Mol_data
    eval(['Mol_data=Mol_total_SecDtime' days(days_i,:) ';']);
    Mol_data_leng=length(Mol_data(:,1));
    figure
    for Mol_data_i=1:Mol_data_leng
        hold on;plot(normalize(Mol_data(Mol_data_i,:),'range'));xlim([0 600])
    end
end

%Check dtime difference during those days
figure
for days_i=1:days_leng
    hold on;
    eval(['plot(SecDtime_matrix_' num2str(days(days_i,:)) '_sum);'])
    xlim([50 1000])
end

figure
for days_i=1:days_leng
    hold on;
    eval(['plot(normalize(SecDtime_matrix_' num2str(days(days_i,:)) '_sum,' char(39) 'range' char(39) '));'])
    xlim([50 1000])
end
