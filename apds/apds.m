%Deal with apd files to use the shape only
clearvars
codefolder=pwd;
solvent='F8T2N2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
srdir=['E:\F8T2400nmCH'];
cd (srdir)


allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);
for len_i=1:1:len
    clearvars -except srdir codefolder solvent len_i len allnames
    datasetname=char(allnames(1,len_i));
    datasetfile=load([srdir '\' datasetname]);
    rowrange=datasetfile.dataset.rowrange;clearvars datasetfile
    disp('Finish load rowrange /n')    
    
    date=regexp(datasetname,'02\d*2020','match');
    file=regexp(datasetname,'\dd\dd\d*','match');
    
    cd([srdir '\apd full'])
    apdfile=dir(['*' date{1} '*' file{1} '.mat']);
    if isempty(apdfile)
        disp(['Wrong apd related to' datasetfile])
    else
        cd(codefolder)
        clearvars apddata apddataresolution
        [apddata,apddataresolution]=PTUim([apdfile.folder '\' apdfile.name]);
    end
    datasource=GetDandABS(apddata,0,'M');
    absolutetime=datasource(:,1);
    dtime=datasource(:,2);
    disp('Finish collect all dtime/n')
    
    
    total_rowrange_leng=0;
    rowrange_leng=length(rowrange);
    for rowrange_i=1:rowrange_leng
       rowrange_rr_leng=length(rowrange(rowrange_i).rr(1,:));
       total_rowrange_leng=rowrange_rr_leng+total_rowrange_leng;
    end
    
    
    SecDtime=cell(total_rowrange_leng,3);flag=1;
    rowrange_leng=length(rowrange);
    for rowrange_i=1:rowrange_leng
       rowrange_rr_leng=length(rowrange(rowrange_i).rr(1,:));
       for rowrange_rr_i=1:rowrange_rr_leng
           SecDtime{flag,1}=dtime(rowrange(rowrange_i).rr(1,rowrange_rr_i):rowrange(rowrange_i).rr(2,rowrange_rr_i),1);
           SecDtime{flag,2}=histcounts(SecDtime{flag,1},1:6252);
           SecDtime{flag,3}=length(SecDtime{flag,1});
           flag=flag+1;
       end
    end
    
cd([srdir '\apd full'])
save(['F8T2 Chloroform 2kDa 400nm ' date{1} ' SecDtime ' file{1} '.mat'],'SecDtime');   
    
end
