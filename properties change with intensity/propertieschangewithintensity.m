%Must have same X% Calculate average wavelength, max wavelength, lifetime
%(choose the first one when calculate in segments), intensity change with
%Intensity.
%E0001 is not settled

clearvars
solvent='F8T2400nmCH';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\F8T2400nmCH'];
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);

apdintensitycombine=[];
lifetimecombine=[];
averagewavelengthcombine=[];
E00wavelengthcombine=[];
Lifindexremove=[];
   
timeave=zeros(99,len);
timemax=zeros(len,99);
timelifetime=zeros(99,len);
%timeE0001=zeros(len,99);
timeintensity=zeros(len,99);
timespectrum=zeros(100,99,len);
timespectrum_normalized=zeros(100,99,len);
place=1;%start to calculate wavelength

for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    try
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %average wavelength change with time%each second,each column
        timeave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with time
        timemax(len_i,:)=datasetfile.dataset.ccdt(maxindex,1);
        %spectrum change with time
        timespectrum(:,:,len_i)=datasetfile.dataset.ccdt(:,3:end);
        timespectrum_normalized(:,:,len_i)=datasetfile.dataset.ccdt(:,3:end)./max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %lifetime change with time
        [newconti_leng,~]=size(datasetfile.dataset.newconti);
        for newconti_i=1:1:newconti_leng
            [~,co_leng]=size(datasetfile.dataset.newconti(newconti_i).co);
            for co_i=1:1:co_leng
                preparemove=datasetfile.dataset.newconti(newconti_i).co(co_i).subco(1,2:end);
                Lifindexremove=cat(2,Lifindexremove,preparemove);
            end
        end
        Lifetime=datasetfile.dataset.scatterplot.lifetime(:,2);
        Lifetime(Lifindexremove,:)=-1;
        timelifetime(:,len_i)=Lifetime;
        %Intensity change with time
        timeintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
%         %E0001 Ratio change with time
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         timeE0001(len_i,:)=E00sum./E01sum;
    catch
        disp([name 'may not have length 99'])
    end
end

edges=430:1:650;
% timeE0001his=zeros(99,20);

intensityedge=min(timeintensity(:)):(max(timeintensity(:))-min(timeintensity(:)))/100:max(timeintensity(:));
intensity_leng=length(intensityedge)-1;
intsp=zeros(100,intensity_leng);intspn=zeros(100,intensity_leng);
timeavehis=zeros(intensity_leng,220);timemaxhis=zeros(intensity_leng,220);timelifetimehis=zeros(intensity_leng,(2500-50)/10);timeintensityhis=zeros(intensity_leng,100);

for intensity_i=1:intensity_leng
    clearvars mol sec 
    [mol,sec]=find(timeintensity>=intensityedge(1,intensity_i) & timeintensity<intensityedge(1,intensity_i+1));
    timeavehis(intensity_i,:)=histcounts(timeave(sec,mol),edges);
    timemaxhis(intensity_i,:)=histcounts(timemax(mol,sec),edges);
    timelifetimehis(intensity_i,:)=histcounts(timelifetime(sec,mol),50:10:2500);
    
    sec_leng=length(sec);
    for sec_i=1:sec_leng
        intsp(:,intensity_i)=intsp(:,intensity_i)+sum(timespectrum(:,sec(sec_i,1),mol(sec_i,1)),2);
        intspn(:,intensity_i)=intspn(:,intensity_i)+normalize(timespectrum(:,sec(sec_i,1),mol(sec_i,1)),1,'range');
    end
end

try
    cd([srdir '/intensity change test/']);
catch
    mkdir([srdir '/intensity change test/']);
    cd([srdir '/intensity change test/']);
end

figure
subplot(1,2,1)
  surf(edges(1,2:end),intensityedge(1,1:end-1),timeavehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['average wavelength change with time ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),intensityedge(1,1:end-1),normalize(timeavehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized average wavelength change with time ' solvent])  
saveas(gcf,[solvent ' average wavelength change with time.jpg']);
  saveas(gcf,[solvent ' average wavelength change with time.fig']);
  disp('Save average wavelength change with time successfully /n');
  close all

figure
subplot(1,2,1)
  surf(edges(1,2:end),intensityedge(1,1:end-1),timemaxhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Maximum Wavelength change with time ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),intensityedge(1,1:end-1),normalize(timemaxhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Maximum wavelength change with time ' solvent])
saveas(gcf,[solvent ' Normalized Maximum change with time.jpg']);
  saveas(gcf,[solvent ' Normalized Maximum change with time.fig']);
  disp('Save Normalized Maximum change with time successfully /n');
  close all
  
figure
subplot(1,2,1)
  surf(60:10:2500,intensityedge(1,1:end-1),timelifetimehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['lifeitme change with time ' solvent])
subplot(1,2,2)
  surf(60:10:2500,intensityedge(1,1:end-1),normalize(timelifetimehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifeitme change with time ' solvent])
saveas(gcf,[solvent ' lifetime change with time.jpg']);
  saveas(gcf,[solvent ' lifetime change with time.fig']);
  disp('Save lifeitme change with time successfully /n');
  close all
    
% figure
% surf(0.05:0.05:1,1:1:99,timeE0001his,'EdgeColor','none');
% colormap(jet)
% view([0 0 1])
%   title(['E0001R change with time ' solvent])
%   saveas(gcf,[solvent ' E0001R change with time.jpg']);
%   saveas(gcf,[solvent ' E0001R change with time.fig']);
%   disp('Save E0001R change with time successfully /n');
%   close all  
  
figure
subplot(1,2,1)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),intsp,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with time ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),normalize(intsp,1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with time ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with time.fig']);
  disp('Save Spectrum (add up) change with time successfully /n');
  close all  
  
figure
subplot(1,2,1)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),intspn,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with time ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),normalize(intspn,1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with time ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with time.fig']);
  disp('Save Spectrum (normalize then add up) change with time successfully /n');
  close all   

 
