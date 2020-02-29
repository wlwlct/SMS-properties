%designed to move dataset, need to put check point to check how good is the
%matching.

% cd('E:\F8T2N2\apd full')
% clearvars
% days=['02012019';'02042019';'02052019'];
% days_leng=length(days(:,1));
% for days_i=1:days_leng
%     day=days(days_i,:);
%     eval(['names_' day '=struct2cell(dir([' char(39) '*' day '*' char(39) ']));']);
%     eval(['current_day=names_' day ';']);
%     oneday_leng=length(current_day(1,:));
%     eval(['name_dtime_' day '=cell(oneday_leng,2);']);
%     for oneday_i=1:oneday_leng
%         file=regexp(current_day(1,oneday_i),'\dd\dd\d*','match');
%         eval(['name_dtime_' day '{oneday_i,1}=file{1,1}{1,1};'])
%         dtime=importdata(current_day{1,oneday_i});
%         eval(['name_dtime_' day '{oneday_i,2}=cell2mat(dtime(:,2));'])
%     end
% end
% 
% dtime_name=who('name_dtime*');
% for i=1:3
%    eval(['current=' dtime_name{i,1} ';'])
%    current_leng=length(current);
%    figure
%    for ii=1:current_leng
%        hold on;clearvars disname
%        eval(['plot(normalize(sum(current{ii,2},1),' char(39) 'range' char(39) '),' char(39) 'DisplayName' char(39) ',' char(39) current{ii,1} char(39) ')'])
%    end
%    xlim([10 500])
%    legend
% end
%%
clearvars
ba=importdata('E:\F8T2N2\apd full\A\F8T2 Chloroform 2kDa N2 02212019 SecDtime 1d2d6.mat');
B=sum(cell2mat(ba(:,2)),1);figure('Position',[-1673 218 560 420]);
date='02222019';
%files={'4d1d10';'4d1d11';'4d1d12';'4d1d13';'4d1d15';'4d1d2';'4d1d4';'4d1d5';'4d1d9'};

move=5;


names=struct2cell(dir(['*' date '*']));
names_leng=length(names(1,:));
files=cell(names_leng,1);
for name_i=1:names_leng
    name=regexp(names{1,name_i},'\dd\dd\d*','match');
    files{name_i,1}=name{1,1};
end
files_leng=length(files(:,1));


for files_i=1:files_leng
    clearvars -except files_i files_leng date files move B ba 
    plot(normalize(B,'range'),'DisplayName','Should be');
    name=dir(['*' date '*' files{files_i,1} '.mat']);
    SecDtime=importdata(name.name);

    SecDtime_leng=length(SecDtime(:,1));SecDtime_ts=cell(SecDtime_leng,3);

    for i=1:SecDtime_leng
        Current_Sec=SecDtime{i,1};
        SecDtime_max=max(Current_Sec);
        SecDtime_ts{i,1}=[Current_Sec(Current_Sec<SecDtime_max-move)+move;Current_Sec(Current_Sec>=SecDtime_max-move)+1-SecDtime_max-move ];
        SecDtime_ts{i,2}=histcounts(SecDtime_ts{i,1},1:6252);
        SecDtime_ts{i,3}=length(SecDtime_ts{i,1}(:,1));
    end

    Sec=sum(cell2mat(SecDtime(:,2)),1);
    Sec_ts=sum(cell2mat(SecDtime_ts(:,2)),1);

    hold on;plot(normalize(Sec,'range'),'DisplayName','original')
    hold on;plot(normalize(Sec_ts,'range'),'DisplayName',['move' num2str(move)]);
    xlim([10 500]);legend;hold off

    
    clearvars -except SecDtime_ts name files_i files_leng date files move B ba 
    SecDtime=SecDtime_ts;
    save([name.name(1:end-4) '_test.mat'],'SecDtime')
end
clearvars;close all;

