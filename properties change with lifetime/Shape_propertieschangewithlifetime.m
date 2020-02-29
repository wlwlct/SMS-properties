%Must have same X% Calculate average wavelength, max wavelength,  intensity change with
%lifetime. %E0001 is not settled

clearvars
solvent='F8T2N2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\F8T2O2'];
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);

Lifindexremove=[];
   
SecDtimeave=zeros(99,len);
SecDtimemax=zeros(len,99);
SecDtimelifetime=zeros(99,len);
%intE0001=zeros(len,99);
SecDtimeintensity=zeros(len,99);
place=22;%start to calculate wavelength
SecDtimespectrum=zeros(100-place+1,99,len);
SecDtimespectrum_normalized=zeros(100-place+1,99,len);

SecDtimeShape=cell(99,len);
SecDtimeShape_bin=zeros(99,len);

edges=450:1:670;
year='2019';

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
        %SecDtimeShape is lifetime shape; SecDtimeShape_bin is the lifetime place
        %reach half the intensity.
        SecDtimeShape_bin(:,len_i)=cellfun(@find50,SecDtime(1:99,2));
        SecDtimeShape(:,len_i)=SecDtime(1:99,2);
        %lifetime change with int
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
        SecDtimelifetime(:,len_i)=lifetime;
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %average wavelength change with int%each second,each column
        SecDtimeave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with int
        SecDtimemax(len_i,:)=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %spectrum change with int
        SecDtimespectrum(:,:,len_i)=datasetfile.dataset.ccdt(place:end,3:end);
        SecDtimespectrum_normalized(:,:,len_i)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range');
        %Intensity change with int
        SecDtimeintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch
        disp([name 'may not have length 99'])
    end
end

%SecDtimeedge=min(SecDtimeShape_bin(:)):(max(SecDtimeShape_bin(:))-min(SecDtimeShape_bin(:)))/100:max(SecDtimeShape_bin(:));
SecDtimeedge=min(SecDtimeShape_bin(:)):max(SecDtimeShape_bin(:));
SecDtimeedge(1,end)=SecDtimeedge(1,end)+1;%include the last element
SecDtimeEdge_leng=length(SecDtimeedge)-1;

SecDtimesp=zeros(100-place+1,SecDtimeEdge_leng);SecDtimespn=zeros(100-place+1,SecDtimeEdge_leng);
SecDtimeavehis=zeros(SecDtimeEdge_leng,220);SecDtimemaxhis=zeros(SecDtimeEdge_leng,220);
SecDtimeintensityhis=zeros(SecDtimeEdge_leng,100);
SecDtimeintensity_edge=min(SecDtimeintensity(:)):(max(SecDtimeintensity(:))-min(SecDtimeintensity(:)))/100:max(SecDtimeintensity(:));
SecDtimelifetime_edge=min(SecDtimelifetime(:)):(max(SecDtimelifetime(:))-min(SecDtimelifetime(:)))/100:max(SecDtimelifetime(:));
SecDtimelifetimehis=zeros(SecDtimeEdge_leng,length(SecDtimelifetime_edge)-1);
SecDtime_sum=zeros(SecDtimeEdge_leng,6251);

for SecDtime_leng_i=1:SecDtimeEdge_leng
    clearvars mol sec 
    [sec,mol]=find((SecDtimeShape_bin>=SecDtimeedge(1,SecDtime_leng_i)) & (SecDtimeShape_bin<SecDtimeedge(1,SecDtime_leng_i+1)));
    sec_leng=length(sec);
    SecDtime_max_prepare=zeros(sec_leng,1); SecDtime_ave_prepare=zeros(sec_leng,1);SecDtime_intensity_prepare=zeros(sec_leng,1);
    SecDtime_lifetime_prepare=zeros(sec_leng,1);
    for sec_i=1:sec_leng
        SecDtimesp(:,SecDtime_leng_i)=SecDtimesp(:,SecDtime_leng_i)+SecDtimespectrum(:,sec(sec_i,1),mol(sec_i,1));
        SecDtimespn(:,SecDtime_leng_i)=SecDtimespn(:,SecDtime_leng_i)+SecDtimespectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
        SecDtime_ave_prepare(sec_i,1)=SecDtimeave(sec(sec_i,1),mol(sec_i,1));
        SecDtime_max_prepare(sec_i,1)=SecDtimemax(mol(sec_i,1),sec(sec_i,1));
        SecDtime_intensity_prepare(sec_i,1)=SecDtimeintensity(mol(sec_i,1),sec(sec_i,1));
        SecDtime_lifetime_prepare(sec_i,1)=SecDtimelifetime(sec(sec_i,1),mol(sec_i,1));
        SecDtime_sum(SecDtime_leng_i,:)=SecDtime_sum(SecDtime_leng_i,:)+cell2mat(SecDtimeShape(sec(sec_i,1),mol(sec_i,1)));
    end
    SecDtimeavehis(SecDtime_leng_i,:)=histcounts(SecDtime_ave_prepare,edges);
    SecDtimemaxhis(SecDtime_leng_i,:)=histcounts(SecDtime_max_prepare,edges);
    SecDtimeintensityhis(SecDtime_leng_i,:)=histcounts(SecDtime_intensity_prepare,SecDtimeintensity_edge);
    SecDtimelifetimehis(SecDtime_leng_i,:)=histcounts(SecDtime_lifetime_prepare,SecDtimelifetime_edge);
end

try
    cd([srdir '/lifetime change/']);
catch
    mkdir([srdir '/lifetime change/']);
    cd([srdir '/lifetime change/']);
end

save([solvent 'SecDtimeShape.mat'],'SecDtimeShape')

figure
subplot(1,2,1)
  surf(edges(1,2:end),SecDtimeedge(1,1:end-1),SecDtimeavehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['average wavelength change with lifetime shape ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),SecDtimeedge(1,1:end-1),normalize(SecDtimeavehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized average wavelength change with lifetime shape ' solvent])  
saveas(gcf,[solvent ' average wavelength change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' average wavelength change with lifetime shape.fig']);
  disp('Save average wavelength change with lifetime shape successfully /n');
  close all

  
figure
subplot(1,2,1)
  surf(edges(1,2:end),SecDtimeedge(1,1:end-1),SecDtimemaxhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Maximum Wavelength change with lifetime shape ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),SecDtimeedge(1,1:end-1),normalize(SecDtimemaxhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Maximum wavelength change with lifetime shape ' solvent])
saveas(gcf,[solvent ' Normalized Maximum change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' Normalized Maximum change with lifetime shape.fig']);
  disp('Save Normalized Maximum change with lifetime shape successfully /n');
  close all
  
figure
subplot(1,2,1)
  surf(SecDtimeintensity_edge(2:end),SecDtimeedge(1,1:end-1),SecDtimeintensityhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Intensity change with lifetime shape ' solvent])
subplot(1,2,2)
  surf(SecDtimeintensity_edge(2:end),SecDtimeedge(1,1:end-1),normalize(SecDtimeintensityhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized intensity change with lifetime shape ' solvent])
saveas(gcf,[solvent ' intensity change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' intenstity change with lifetime shape.fig']);
  disp('Save intensity change with lifetime shape successfully /n');
  close all
    
  
figure
subplot(1,2,1)
  surf(SecDtimelifetime_edge(2:end),SecDtimeedge(1,1:end-1),SecDtimelifetimehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['lifetime number change with lifetime shape ' solvent])
subplot(1,2,2)
  surf(SecDtimelifetime_edge(2:end),SecDtimeedge(1,1:end-1),normalize(SecDtimelifetimehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifetime number change with lifetime shape ' solvent])
saveas(gcf,[solvent ' lifetime number change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' lifetime number change with lifetime shape.fig']);
  disp('Save lifetime number change with lifetime shape successfully /n');
  close all
  
figure
subplot(2,2,1)
  histogram(SecDtimeShape_bin(:),min(SecDtimeShape_bin(:)):max(SecDtimeShape_bin(:)))
  title(['Distribution of lifetime shape ' solvent])
subplot(2,2,2)
  surf(8*(50:1000)/1000,SecDtimeedge(1,1:end),[SecDtime_sum(:,50:1000);zeros(1,951)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['lifetime number change with lifetime shape ' solvent])
subplot(2,2,3)
  surf(8*(50:1000)/1000,SecDtimeedge(1,1:end),normalize([SecDtime_sum(:,50:1000); zeros(1,951)],2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifetime dtime shape change with lifetime shape' solvent])
subplot(2,2,4)
  contourf(8*(50:1000)/1000,SecDtimeedge(1,1:end),normalize([SecDtime_sum(:,50:1000); zeros(1,951)],2,'range'),1);colormap(jet);view([0 0 1]);
  title(['Normalized lifetime dtime shape change with lifetime shape' solvent]) 
saveas(gcf,[solvent ' lifetime dtime shape change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' lifetime dtime shape change with lifetime shape.fig']);
  disp('Save lifetime dtime shape change with lifetime shape successfully /n');
  close all
  
  
% figure
% surf(0.05:0.05:1,1:1:99,intE0001his,'EdgeColor','none');
% colormap(jet)
% view([0 0 1])
%   title(['E0001R change with int ' solvent])
%   saveas(gcf,[solvent ' E0001R change with int.jpg']);
%   saveas(gcf,[solvent ' E0001R change with int.fig']);
%   disp('Save E0001R change with int successfully /n');
%   close all  
  
figure
subplot(1,2,1)
  surf(SecDtimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),[SecDtimesp(:,:) zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with lifetime shape ' solvent])
subplot(1,2,2)
  surf(SecDtimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),normalize([SecDtimesp(:,:) zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with lifetime shape ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with lifetime shape.fig']);
  disp('Save Spectrum (add up) change with lifetime shape successfully /n');
  close all  
  
figure
subplot(1,2,1)
  surf(SecDtimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),[SecDtimespn(:,:) zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with lifetime shape ' solvent])
subplot(1,2,2)
  surf(SecDtimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),normalize([SecDtimespn(:,:) zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with lifetime shape ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with lifetime shape.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with lifetime shape.fig']);
  disp('Save Spectrum (normalize then add up) change with lifetime shape successfully /n');
  close all   
 
  
  function bin=find50(rowshape)
    rowshape=rowshape(50:2000);
    intensity=sum(rowshape,2);
    norm_rowshape=rowshape./intensity;
    rowshape_leng=length(rowshape);
    halfplace=ones(1,ceil(rowshape_leng/2))*(-1);
    for i=1:ceil(rowshape_leng/2)
        halfplace(1,i)=sum(norm_rowshape(1:i),2);
    end
    [~,bin]=min(abs(halfplace-0.5));
  end
