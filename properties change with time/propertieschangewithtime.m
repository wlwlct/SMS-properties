%Must have same X% Calculate average wavelength, max wavelength, lifetime
%(choose the first one when calculate in segments), intensity change with
%time.
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
   
timeave=zeros(99,len);
timemax=zeros(len,99);
timelifetime=zeros(99,len);
timeE0001=zeros(len,99);
timeintensity=zeros(len,99);
timespectrum=zeros(100,99);  
timespectrum_normalized=zeros(100,99);   
place=1;%start to calculate wavelength

for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    try
        
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        wavelengthindex=datasetfile.dataset.ccdt(maxindex,1);
        %average wavelength change with time%each second,each column
        timeave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with time
        timemax(len_i,:)=wavelengthindex;
        %spectrum change with time
        timespectrum=datasetfile.dataset.ccdt(:,3:end)+timespectrum;
        timespectrum_normalized=datasetfile.dataset.ccdt(:,3:end)./max(datasetfile.dataset.ccdt(place:end,3:end),[],1)+timespectrum_normalized;
        %lifetime change with time
        Lifindexremove=[];
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
%%
edges=430:1:650;intensityedge=min(timeintensity(:)):(max(timeintensity(:))-min(timeintensity(:)))/100:max(timeintensity(:));
int10000=1:100:10000;
timeavehis=zeros(99,220);timemaxhis=zeros(99,220);timelifetimehis=zeros(99,(2500-50)/10);timeintensityhis=zeros(99,100);
% timeE0001his=zeros(99,20);

%histogramize
for i=1:1:99
    timeavehis(i,:)=histcounts(timeave(i,:),edges);
    timemaxhis(i,:)=histcounts(timemax(:,i),edges);
    timelifetimehis(i,:)=histcounts(timelifetime(i,:),50:10:2500);
    timeintensityhis(i,:)=histcounts(timeintensity(:,i),intensityedge);
    timeint10000his(i,:)=histcounts(timeintensity(:,i),int10000);
%     timeE0001his(i,:)=histcounts(timeE0001(:,i),0:0.05:1);
end

meanintensity=mean(timeintensity,1);
stdintensity=std(timeintensity,1);

try
    cd([srdir '/time change/']);
catch
    mkdir([srdir '/time change/']);
    cd([srdir '/time change/']);
end

figure
subplot(1,2,1)
  surf(edges(1,2:end),1:1:99,timeavehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Average wavelength change with time ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),1:1:99,normalize(timeavehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized average wavelength change with time ' solvent])  
saveas(gcf,[solvent ' average wavelength change with time.jpg']);
  saveas(gcf,[solvent ' average wavelength change with time.fig']);
  disp('Save average wavelength change with time successfully /n');
  close all

figure
subplot(1,2,1)
  surf(edges(1,2:end),1:1:99,timemaxhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Maximum Wavelength change with time ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),1:1:99,normalize(timemaxhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized maximum wavelength change with time ' solvent])
saveas(gcf,[solvent ' Normalized Maximum change with time.jpg']);
  saveas(gcf,[solvent ' Normalized Maximum change with time.fig']);
  disp('Save Normalized Maximum change with time successfully /n');
  close all
  
figure
subplot(1,2,1)
  surf(60:10:2500,1:1:99,timelifetimehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Lifeitme change with time ' solvent])
subplot(1,2,2)
  surf(60:10:2500,1:1:99,normalize(timelifetimehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifeitme change with time ' solvent])
saveas(gcf,[solvent ' lifetime change with time.jpg']);
  saveas(gcf,[solvent ' lifetime change with time.fig']);
  disp('Save lifeitme change with time successfully /n');
  close all
  
figure
subplot(1,2,1)
  surf(intensityedge(1,2:end),1:1:99,timeintensityhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Intensity change with time ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,2:end),1:1:99,normalize(timeintensityhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized intensity change with time ' solvent])
saveas(gcf,[solvent ' Intensity change with time.jpg']);
  saveas(gcf,[solvent ' Intensity change with time.fig']);
  disp('Save Intensity change with time successfully /n');
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
  surf(1:1:99,datasetfile.dataset.ccdt(:,1),timespectrum,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with time ' solvent])
subplot(1,2,2)
  surf(1:1:99,datasetfile.dataset.ccdt(:,1),normalize(timespectrum,1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with time ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with time.fig']);
  disp('Save Spectrum (add up) change with time successfully /n');
  close all  
  
figure
subplot(1,2,1)
  surf(1:1:99,datasetfile.dataset.ccdt(:,1),timespectrum_normalized,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with time ' solvent])
subplot(1,2,2)
  surf(1:1:99,datasetfile.dataset.ccdt(:,1),normalize(timespectrum_normalized,1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with time ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with time.fig']);
  disp('Save Spectrum (normalize then add up) change with time successfully /n');
  close all   

 
figure
subplot(1,2,1)
  surf(int10000(1,2:end),1:1:99,timeint10000his,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Intensity 10000 change with time ' solvent])
subplot(1,2,2)
  surf(int10000(1,2:end),1:1:99,normalize(timeint10000his,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized intensity 10000 change with time ' solvent])
  saveas(gcf,[solvent ' Intensity 10000 change with time.jpg']);
  saveas(gcf,[solvent ' Intensity 10000 change with time.fig']);
  disp('Save Intensity 10000 change with time successfully /n');
  close all  
  

 figure;yyaxis left; plot(meanintensity);ylabel('mean');yyaxis right;plot(stdintensity); ylabel('standard deviation')
  xlabel('time');title([solvent 'Intensity change with time'])
  saveas(gcf,[solvent ' Intensity plot change with time.jpg']);
  saveas(gcf,[solvent ' Intensity plot change with time.fig']);
  disp('Save Intensity plot change with time successfully /n');
  close all  
