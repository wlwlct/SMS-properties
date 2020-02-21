%Must have same X% Calculate average wavelength, max wavelength, lifetime
%(choose the first one when calculate in segments), intensity change with
%Intensity. Overall avewav, maxave, lifetime, spectrum(normalized and not normalized.)
%E0001 is not settled

clearvars
solvent='F8T2O2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\F8T2400nmCH'];
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);

Lifindexremove=[];
   
intave=zeros(99,len);
intmax=zeros(len,99);
intlifetime=zeros(99,len);
%intE0001=zeros(len,99);
intintensity=zeros(len,99);
intspectrum=zeros(100,99,len);
intspectrum_normalized=zeros(100,99,len);
place=22;%start to calculate wavelength
edges=450:1:670;

for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    try
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %average wavelength change with int%each second,each column
        intave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with int
        intmax(len_i,:)=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %spectrum change with int
        intspectrum(:,:,len_i)=datasetfile.dataset.ccdt(:,3:end);
        intspectrum_normalized(:,:,len_i)=datasetfile.dataset.ccdt(:,3:end)./max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
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
        intlifetime(:,len_i)=lifetime;
        %Intensity change with int
        intintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch
        disp([name 'may not have length 99'])
    end
end

%edges=430:1:650;
% intE0001his=zeros(99,20);

intensityedge=min(intintensity(:)):(max(intintensity(:))-min(intintensity(:)))/100:max(intintensity(:));
intensity_leng=length(intensityedge)-1;
intsp=zeros(100,intensity_leng);intspn=zeros(100,intensity_leng);
intavehis=zeros(intensity_leng,220);intmaxhis=zeros(intensity_leng,220);intlifetimehis=zeros(intensity_leng,(2500-50)/10);intintensityhis=zeros(intensity_leng,100);

for intensity_i=1:intensity_leng
    clearvars mol sec 
    [mol,sec]=find((intintensity>=intensityedge(1,intensity_i)) & (intintensity<intensityedge(1,intensity_i+1)));
    sec_leng=length(sec);
    intmax_prepare=zeros(sec_leng,1); intave_prepare=zeros(sec_leng,1);intlifetime_prepare=zeros(sec_leng,1);
    for sec_i=1:sec_leng
        intsp(:,intensity_i)=intsp(:,intensity_i)+intspectrum(:,sec(sec_i,1),mol(sec_i,1));
        intspn(:,intensity_i)=intspn(:,intensity_i)+intspectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
        intave_prepare(sec_i,1)=intave(sec(sec_i,1),mol(sec_i,1));
        intmax_prepare(sec_i,1)=intmax(mol(sec_i,1),sec(sec_i,1));
        intlifetime_prepare(sec_i,1)=intlifetime(sec(sec_i,1),mol(sec_i,1));
    end
    intavehis(intensity_i,:)=histcounts(intave_prepare,edges);
    intmaxhis(intensity_i,:)=histcounts(intmax_prepare,edges);
    intlifetimehis(intensity_i,:)=histcounts(intlifetime_prepare,50:10:2500);
end

try
    cd([srdir '/intensity change/']);
catch
    mkdir([srdir '/intensity change/']);
    cd([srdir '/intensity change/']);
end

figure
subplot(1,2,1)
  surf(edges(1,2:end),intensityedge(1,1:end-1),intavehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['average wavelength change with int ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),intensityedge(1,1:end-1),normalize(intavehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized average wavelength change with int ' solvent])  
saveas(gcf,[solvent ' average wavelength change with int.jpg']);
  saveas(gcf,[solvent ' average wavelength change with int.fig']);
  disp('Save average wavelength change with int successfully /n');
  close all

  figure
subplot(2,2,1)
  plot(edges(1,2:end),sum(intmaxhis,1),'LineWidth',3);
  title(['Overall max wavelength distribution ' solvent]) 
  hold on;plot(edges(1,2:end),sum(intavehis,1),'LineWidth',3);
  title(['Overall average wavelength distribution ' solvent])
subplot(2,2,2)
  histogram(intintensity,intensityedge);
  title(['Overall intensity distribution ' solvent])
subplot(2,2,3)
  plot(60:10:2500,sum(intlifetimehis,1),'LineWidth',3);
  title(['Overall lifetime distribution ' solvent]) 
subplot(2,2,4)
  yyaxis left; plot(datasetfile.dataset.ccdt(:,1),sum(intsp,2),'LineWidth',3);
  title(['Overall spectrum add up ' solvent])
  yyaxis right;plot(datasetfile.dataset.ccdt(:,1),sum(intspn,2),'LineWidth',3);
  title(['Overall normalized spectrum add up ' solvent])
saveas(gcf,[solvent ' all add up.jpg']);
  saveas(gcf,[solvent ' all add up.fig']);
  disp('Save all add up successfully /n');
  close all

  
figure
subplot(1,2,1)
  surf(edges(1,2:end),intensityedge(1,1:end-1),intmaxhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Maximum Wavelength change with int ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),intensityedge(1,1:end-1),normalize(intmaxhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Maximum wavelength change with int ' solvent])
saveas(gcf,[solvent ' Normalized Maximum change with int.jpg']);
  saveas(gcf,[solvent ' Normalized Maximum change with int.fig']);
  disp('Save Normalized Maximum change with int successfully /n');
  close all
  
figure
subplot(1,2,1)
  surf(60:10:2500,intensityedge(1,1:end-1),intlifetimehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['lifeitme change with int ' solvent])
subplot(1,2,2)
  surf(60:10:2500,intensityedge(1,1:end-1),normalize(intlifetimehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifeitme change with int ' solvent])
saveas(gcf,[solvent ' lifetime change with int.jpg']);
  saveas(gcf,[solvent ' lifetime change with int.fig']);
  disp('Save lifeitme change with int successfully /n');
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
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),intsp,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with int ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),intsp/max(intsp(place:end,:),[],1),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with int ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with int.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with int.fig']);
  disp('Save Spectrum (add up) change with int successfully /n');
  close all  
  
figure
subplot(1,2,1)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),intspn,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with int ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end-1),datasetfile.dataset.ccdt(:,1),intspn/max(intspn(place:end,:),[],1),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with int ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with int.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with int.fig']);
  disp('Save Spectrum (normalize then add up) change with int successfully /n');
  close all   
 
