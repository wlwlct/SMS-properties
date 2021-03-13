%Parameter setting
clearvars;solvent='F8T2N2';
codefolder=pwd;
%srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\F8Se2 July\' 'F8Se2O2'];
srdir=['E:\F8T2400nmCH'];
cd (srdir)
allnames=struct2cell(dir( '*dataset*.mat'));
[~,len]=size(allnames);

%Find ordered by max; ordered by ratio in three different range.
edges=450:1:670;
year='2020';
place=1;%start to calculate wavelength
max_secs=200;%choose according to the max number of seconds in each files.

Lifindexremove=[];
spectralifetime=zeros(max_secs,len);
spectraintensity=zeros(len,max_secs);
spectraspectrum=zeros(100-place+1,max_secs,len);
spectraspectrum_normalized=zeros(100-place+1,max_secs,len);
SecDtimeintensity=cell(max_secs,len);
spectraspectrum_diff=zeros(100-place+1,max_secs-1,len);
spectraspectrum_diff_std=zeros(len,max_secs-1);
spectra_stage=zeros(100-place+1,max_secs,len);
spectra_stage_ratio=zeros(max_secs,len);
%import data; only contain raw data;
for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    
    date=regexp(name,['\d*' year],'match');
    file=regexp(name,'\dd\dd\d*','match');
    cd([srdir '/apd full'])
    Secfile=dir(['*' date{1} '*SecDtime*' file{1} '.mat']);
    SecDtime=importdata(Secfile.name);
    
    try
        %max wavelength change with spectra
        %spectrum change with spectra
        spectraspectrum_leng=length(datasetfile.dataset.ccdt(1,3:end));
        spectraspectrum(:,1:spectraspectrum_leng,len_i)=datasetfile.dataset.ccdt(place:end,3:end);
        for ii=1:max_secs
            try
                clearvars A
                cd(codefolder)
                [~,A.eff_fit,~,A.numst,~]=Traceson(spectraspectrum(:,ii,len_i),codefolder);
                if ~isempty(A)
                    if length(A.eff_fit(:,1))<A.numst;efffit=A.eff_fit(1,:);else;efffit=A.eff_fit(A.numst,:);end
                    spectra_stage(:,ii,len_i)=transpose(efffit);
                    spectra_stage_ratio(ii,len_i)=max(spectra_stage(:,ii,len_i))/min(spectra_stage(:,ii,len_i));
                end
            catch
            end
        end
        spectraspectrum_normalized(:,1:spectraspectrum_leng,len_i)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range');
        %lifetime change with spectra
        [newconti_leng,~]=size(datasetfile.dataset.newconti);
        for newconti_i=1:1:newconti_leng
            [~,co_leng]=size(datasetfile.dataset.newconti(newconti_i).co);
            for co_i=1:1:co_leng
                preparemove=datasetfile.dataset.newconti(newconti_i).co(co_i).subco(1,2:end);
                Lifindexremove=cat(2,Lifindexremove,preparemove);
            end
        end
        lifetime=datasetfile.dataset.scatterplot.lifetime(:,2);
        lifetime(Lifindexremove,:)=-1;
        spectralifetime_leng=length(lifetime);
        spectralifetime(1:spectralifetime_leng,len_i)=lifetime;
        %Intensity change with spectra
        spectraintensity_leng=length(datasetfile.dataset.scatterplot.intensity(1,:));
        spectraintensity(len_i,1:spectraintensity_leng)=datasetfile.dataset.scatterplot.intensity(1,:);
        %dtime change with spectra
        SecDtimeintensity_leng=length(SecDtime(:,2));
        SecDtimeintensity(1:SecDtimeintensity_leng,len_i)=SecDtime(:,2);
        %differentce in spectra;jitter
        spectraspectrum_diff_leng=length(datasetfile.dataset.ccdt(1,3:end))-1;
        spectraspectrum_diff(:,1:spectraspectrum_diff_leng,len_i)=diff(datasetfile.dataset.ccdt(place:end,3:end),1,2);
        spectraspectrum_diff_std(len_i,1:spectraspectrum_diff_leng)=std(spectraspectrum_diff(:,1:spectraspectrum_diff_leng,len_i),0,1);
    catch ME
        rethrow(ME)
        disp([name 'may not have length 99'])
    end
end

%Remove some of the spectrum that has no shape
% spectraspectrum% spectra_stage
% spectraspectrum_diff% spectraspectrum_diff_std
% spectra_stage_ratio% spectraspectrum_normalized
% spectralifetime% spectraintensity% SecDtimeintensity

%spectramax before sign big value to bad spectra
spectramax=zeros(len,max_secs);
for len_i=1:len
    [~,spectramax_loc]=max(spectraspectrum(:,:,len_i),[],1);
    spectramax(len_i,:)=transpose(datasetfile.dataset.ccdt(spectramax_loc+place-1,1));
end

for ra=[1.6 1.7 1.75]
[Bad_sec,Bad_mol]=find(spectra_stage_ratio<ra);
Bad_leng=length(Bad_sec);max_int=max(spectraspectrum(:));


for i=1:Bad_leng
    spectraspectrum(end-5:end,Bad_sec(i,1),Bad_mol(i,1))=100+max_int;
    spectraspectrum_normalized(end-5:end,Bad_sec(i,1),Bad_mol(i,1))=1.1;
end

%spectramax before sign big value to bad spectra
spectramax_smooth=zeros(len,max_secs);
for len_i=1:len
    clearvars spectramax_smooth_loc
    [~,spectramax_smooth_loc]=max(smoothdata(spectraspectrum(:,:,len_i),1,'gaussian',8),[],1);
    spectramax_smooth(len_i,:)=transpose(datasetfile.dataset.ccdt(spectramax_smooth_loc+place-1,1));
end

%cut anything related to the last second, it would not related to spectrum
%change or jitter.

spectramax_smooth=spectramax_smooth(:,1:max_secs-1);

%Sort spectra based on maxwavelength. Inside each max wavelength, the
%spectra is sorted or not based on the similarity.
spectrum_edge_leng=length(edges)-1;
spectra_max_prepare=cell(spectrum_edge_leng,1);
spectra_next_prepare=cell(spectrum_edge_leng,1);
spectra_intensity_prepare=cell(spectrum_edge_leng,1);
spectra_Dtime_prepare=cell(spectrum_edge_leng,1);

%Try not to sort in the first place
for spectrum_edge_i=1:spectrum_edge_leng 
    clearvars mol sec
    [mol,sec]=find((spectramax_smooth>=edges(1,spectrum_edge_i)) & (spectramax_smooth<edges(1,spectrum_edge_i+1)));
    sec_leng=length(sec); 
        spectra_max_prepare{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
        spectra_next_prepare{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
        spectra_intensity_prepare{spectrum_edge_i,1}=zeros(1,sec_leng);
        spectra_Dtime_prepare{spectrum_edge_i,1}=cell(1,sec_leng);
    for sec_i=1:sec_leng
        spectra_max_prepare{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
        spectra_next_prepare{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized(:,sec(sec_i,1)+1,mol(sec_i,1));
        spectra_intensity_prepare{spectrum_edge_i}(1,sec_i)=spectraintensity(mol(sec_i,1),sec(sec_i,1));
        spectra_Dtime_prepare{spectrum_edge_i}{1,sec_i}=SecDtimeintensity{sec(sec_i,1),mol(sec_i,1)};
    end
end

try
    cd([srdir '/spectra change/']);
catch
    mkdir([srdir '/spectra change/']);
    cd([srdir '/spectra change/']);
end

%average spectra shape in each range,mesh current next and difference...
spectra_max_average=zeros(100-place+1,220);
spectra_next_average=zeros(100-place+1,220);
spectra_intensity_average=zeros(2,220);%mean and std
spectra_Dtime_average=zeros(6251,220);
for i=1:220
    if ~isempty(spectra_max_prepare{i,1})
        spectra_max_average(:,i)=mean(spectra_max_prepare{i,1},2);
        spectra_next_average(:,i)=mean(spectra_next_prepare{i,1},2);
        spectra_intensity_average(1,i)=mean(spectra_intensity_prepare{i,1},2);
        spectra_intensity_average(2,i)=std(spectra_intensity_prepare{i,1},0,2);
        
        clearvars Dtime_leng Dtime_matrix
        Dtime_leng=length(spectra_Dtime_prepare{i,1}(1,:));
        Dtime_matrix=zeros(Dtime_leng,6251);
        for Dtime_i=1:Dtime_leng
            Dtime_matrix(Dtime_i,:)=normalize(spectra_Dtime_prepare{i,1}{1,Dtime_i},2,'range');
        end
        spectra_Dtime_average(:,i)=transpose(mean(Dtime_matrix,1));
    end
end

if ra<0
    ra=complex(0,ra);
end

%spectrum distribution stuff
figure
hold on; histogram(spectramax(:),400:900)
hold on; histogram(spectramax_smooth(:),400:900)
legend('spectra_m_a_x','spectra_m_a_x_s_o_o_t_h')
title(ra)
saveas(gcf,[num2str(ra) ' dis of max and smooth max.jpg']);
saveas(gcf,[num2str(ra) ' dis of max and smooth max.fig']);
close all

%plot average max spectrum; next spectrum; intensity; Dtime;
figure('Position',[0,0,729,554]);
subplot(2,2,1);mesh(edges,datasetfile.dataset.ccdt(place:end,1),normalize([spectra_max_average zeros(100-place+1,1)],1,'range'));
    view([0 0 1]); colormap(jet);title('current');ylabel('Wavelength (nm)');xlabel('Max Wavelength (nm)');
subplot(2,2,2);mesh(edges,datasetfile.dataset.ccdt(place:end,1),normalize([spectra_next_average zeros(100-place+1,1)],1,'range'));
    view([0 0 1]); colormap(jet);title('next');ylabel('Wavelength (nm)');xlabel('Max Wavelength (nm)');
subplot(2,2,3);yyaxis left;plot(edges(1:220),spectra_intensity_average(1,:));ylabel('Mean Intensity')
    yyaxis right;plot(edges(1:220),spectra_intensity_average(2,:));ylabel('Std Intensity');xlabel('Max Wavelength (nm)')
subplot(2,2,4);mesh(edges,(1:6251)*8/1000,[normalize(spectra_Dtime_average,1,'range') zeros(6251,1)]);
    view([0 0 1]); colormap(jet);title('Dtime');ylim([1 5]);xlabel('Max Wavelength (nm)');ylabel('Dtime (ns)');

saveas(gcf,[num2str(ra) solvent ' Dtime int with spectra with blank.fig']);
saveas(gcf,[num2str(ra) ' Dtime int with spectra with blank.jpg']);
close all

loc=find(any(spectra_max_average));loc_leng=length(loc(1,:));loc_name=cellfun(@num2str,num2cell(edges(1,loc)),'UniformOutput',false);
figure('Position',[0,0,729,554]);
subplot(2,2,1);mesh(1:loc_leng+1,datasetfile.dataset.ccdt(place:end,1),normalize([spectra_max_average(:,loc) zeros(100-place+1,1)],1,'range'));
    xticks(1:loc_leng+1);xticklabels([loc_name '0']);
    view([0 0 1]); colormap(jet);title('current');ylabel('Wavelength (nm)');xlabel('Max Wavelength (nm)');
subplot(2,2,2);mesh(1:loc_leng+1,datasetfile.dataset.ccdt(place:end,1),normalize([spectra_next_average(:,loc) zeros(100-place+1,1)],1,'range'));
    xticks(1:loc_leng+1);xticklabels([loc_name '0']);
    view([0 0 1]); colormap(jet);title('next');ylabel('Wavelength (nm)');xlabel('Max Wavelength (nm)');
subplot(2,2,3);yyaxis left;plot(1:loc_leng,spectra_intensity_average(1,loc));ylabel('Mean Intensity')
    yyaxis right;plot(1:loc_leng,spectra_intensity_average(2,loc));ylabel('Std Intensity');xlabel('Max Wavelength (nm)')
    xticks(1:loc_leng);xticklabels(loc_name);
subplot(2,2,4);mesh(1:loc_leng+1,(1:6251)*8/1000,[normalize(spectra_Dtime_average(:,loc),1,'range') zeros(6251,1)]);
    view([0 0 1]); colormap(jet);title('Dtime');ylim([1 5]);xlabel('Max Wavelength (nm)');ylabel('Dtime (ns)');
    xticks(1:loc_leng+1);xticklabels([loc_name '0']);

saveas(gcf,[num2str(ra) solvent ' Dtime int with spectra.fig']);
saveas(gcf,[num2str(ra) solvent ' Dtime int with spectra.jpg']);
close all
end