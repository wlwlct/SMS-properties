%Parameter setting
clearvars;solvent='F8T2SMS400';
codefolder=pwd;
%srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
srdir=['E:\F8T2400nmCH'];
cd (srdir)
allnames=struct2cell(dir( '*dataset*.mat'));
[~,len]=size(allnames);

%Find ordered by max; ordered by ratio in three different range.
edges=400:1:670;
year='2020';
place=1;%start to calculate wavelength
max_secs=200;%choose according to the max number of seconds in each files.


Lifindexremove=[];
spectramax=zeros(len,max_secs);
spectraintensity=zeros(len,max_secs);

spectraspectrum=zeros(100-place+1,max_secs,len);
spectraspectrum_normalized=zeros(100-place+1,max_secs,len);
SecDtimeintensity=cell(max_secs,len);
spectraspectrum_diff=zeros(100-place+1,max_secs,len);
spectraspectrum_diff_std=zeros(len,max_secs);
spectra_stage=zeros(100-place+1,max_secs,len);
spectra_stage_ratio=zeros(max_secs,len);

%import data; only contain raw data;
for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp(['Finish load file /n',num2str(len_i)])
    
    date=regexp(name,['\d*' year],'match');
    file=regexp(name,'\dd\dd\d*','match');
    cd([srdir '/apd full'])
    Secfile=dir(['*' date{1} '*SecDtime*' file{1} '.mat']);
    SecDtime=importdata(Secfile.name);
    
    try
        [~,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %max wavelength change with spectra
        spectramax_leng=length(datasetfile.dataset.ccdt(maxindex+place-1,1));
        spectramax(len_i,1:spectramax_leng)=datasetfile.dataset.ccdt(maxindex+place-1,1);
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
        %Intensity change with spectra
        spectraintensity_leng=length(datasetfile.dataset.scatterplot.intensity(1,:));
        spectraintensity(len_i,1:spectraintensity_leng)=datasetfile.dataset.scatterplot.intensity(1,:);
        %dtime change with spectra
        %SecDtimeintensity(:,len_i)=SecDtime(1:max_secs,2);
        %differentce in spectra;jitter
        spectraspectrum_diff_leng=length(datasetfile.dataset.ccdt(1,3:end))-1;
        spectraspectrum_diff(:,1:spectraspectrum_diff_leng,len_i)=diff(datasetfile.dataset.ccdt(place:end,3:end),1,2);
        spectraspectrum_diff_std(len_i,:)=std(spectraspectrum_diff(:,:,len_i),0,1);
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch ME
        rethrow(ME)
        disp([name 'may not have length max_secs'])
    end
end

%Remove some of the spectrum that has no shape
% spectraspectrum% spectra_stage
% spectraspectrum_diff% spectraspectrum_diff_std
% spectra_stage_ratio% spectraspectrum_normalized% spectramax
% spectraintensity% SecDtimeintensity

[Bad_sec,Bad_mol]=find(spectra_stage_ratio<-1);
Bad_leng=length(Bad_sec);max_int=max(spectraspectrum(:));
for i=1:Bad_leng
    spectraspectrum(end-5:end,Bad_sec(i,1),Bad_mol(i,1))=100+max_int;
    spectraspectrum_normalized(end-5:end,Bad_sec(i,1),Bad_mol(i,1))=1.1;
end
spectramax_smooth=zeros(len,max_secs);
for len_i=1:len
    clearvars spectramax_smooth_loc
[~,spectramax_smooth_loc]=max(smoothdata(spectraspectrum(:,:,len_i),1,'gaussian',8),[],1);
spectramax_smooth(len_i,:)=transpose(datasetfile.dataset.ccdt(spectramax_smooth_loc,1));
end

%cut anything related to the last second, it would not related to spectrum
%change or jitter.
spectramax=spectramax(:,1:max_secs);
spectramax_smooth=spectramax_smooth(:,1:max_secs);
spectraspectrum_normalized_smooth=smoothdata(spectraspectrum_normalized,1,'gaussian',8);

%Sort spectra based on maxwavelength. Inside each max wavelength, the
%spectra is sorted or not based on the similarity.
spectrum_edge_leng=length(edges)-1;
spectra_max_prepare=cell(spectrum_edge_leng,1);spectra_max_prepare_smooth=cell(spectrum_edge_leng,1);
spectra_next_intensity_prepare=cell(spectrum_edge_leng,1);
spectra_jitter_prepare=cell(spectrum_edge_leng,1);
spectra_jitter_prepare_normalize=cell(spectrum_edge_leng,1);
spectra_next_prepare=cell(spectrum_edge_leng,1);spectra_next_prepare_smooth=cell(spectrum_edge_leng,1);
spectra_difference_prepare=cell(spectrum_edge_leng,1);

%Try not to sort in the first place
for spectrum_edge_i=1:spectrum_edge_leng 
    clearvars mol sec
    [mol,sec]=find((spectramax_smooth>=edges(1,spectrum_edge_i)) & (spectramax_smooth<edges(1,spectrum_edge_i+1)));
    sec_leng=length(sec); 
        spectra_max_prepare{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
        spectra_max_prepare_smooth{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
        spectra_next_prepare{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
        spectra_next_prepare_smooth{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
        spectra_difference_prepare{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
    for sec_i=1:sec_leng
        spectra_max_prepare{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
        spectra_max_prepare_smooth{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized_smooth(:,sec(sec_i,1),mol(sec_i,1));
%         spectra_next_intensity_prepare{spectrum_edge_i}(1,sec_i)=spectraintensity(mol(sec_i,1),sec(sec_i,1)+1);
        spectra_jitter_prepare{spectrum_edge_i}(1,sec_i)=spectraspectrum_diff_std(mol(sec_i,1),sec(sec_i,1));
%        spectra_jitter_prepare_normalize{spectrum_edge_i}(1,sec_i)=spectraspectrum_std_normalize(mol(sec_i,1),sec(sec_i,1));
    end
end
wl=datasetfile.dataset.ccdt(place:end,1);
%%
try
    cd([srdir '/spectra change/']);
catch
    mkdir([srdir '/spectra change/']);
    cd([srdir '/spectra change/']);
end
%%
clearvars -except wl spectraspectrum_normalized*
[~,sp]=max(spectraspectrum_normalized_smooth,[],1);
sp=squeeze(sp);%verticle spectrum; horizontal each molecule
sp=wl(sp);
testL530=double((sp>=520) & (sp<700));
testS530=double(sp<520);
testL700=double(sp>=700);
X=repmat(1:94,200,1);
Y=transpose(repmat(1:200,94,1));

wid=50;
figure('Position',[349,545,1365,321])
hold on;scatter3(X(:),Y(:),testL700(:),wid,'Marker','square','MarkerFaceColor',[255, 234, 167]./256,'MarkerEdgeColor','none')
hold on;scatter3(X(:),Y(:),testS530(:),wid,'Marker','square','MarkerFaceColor',[9, 132, 227]./256,'MarkerEdgeColor','none')
hold on;scatter3(X(:),Y(:),testL530(:),wid,'Marker','square','MarkerFaceColor',[0, 184, 148]./256,'MarkerEdgeColor','none')
zlim([0.5 1.5]);view([0 0 1])

box on
ax=gca;
ax.LineWidth=1.5;
ax.FontSize=12;
ax.FontWeight='bold'
xlim([0 113])
ylim([0 100])
ylabel('Experimental Time (s)')
xlabel('ID of Each Molecule')
legend('Dim State','E_m_a_x<520 nm','E_m_a_x\geq520 nm')
