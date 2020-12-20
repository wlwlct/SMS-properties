%Must have same X% Calculate average wavelength, max wavelength, lifetime
%(choose the first one when calculate in segments), intensity change with
%time.
%E0001 is not settled

clearvars
solvent='F8Se2O2';
%srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
srdir=['E:\F8Se2 July\' solvent];
%srdir=['E:\F8T2400nmCH'];
max_secs=200;%choose according to the max number of seconds in each files.
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);

apdintensitycombine=[];
lifetimecombine=[];
averagewavelengthcombine=[];
E00wavelengthcombine=[];
   
timeave=zeros(max_secs,len);
timemax=zeros(len,max_secs);
timelifetime=zeros(max_secs,len);
timeE0001=zeros(len,max_secs);
timeintensity=zeros(len,max_secs);
place=22;%start to calculate wavelength
timespectrum=zeros(100-place+1,max_secs);  
timespectrum_normalized=zeros(100-place+1,max_secs);   
SecDtimeintensity=cell(max_secs,len);
edges=450:1:670;
year='2020';

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
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        wavelengthindex=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %average wavelength change with time%each second,each column
        timeave_leng=length(datasetfile.dataset.scatterplot.spectrum(:,1));
        timeave(1:timeave_leng,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with time
        timemax_leng=length(wavelengthindex(:));
        timemax(len_i,1:timemax_leng)=wavelengthindex(:);
        %spectrum change with time
        timespectrum_leng=length(datasetfile.dataset.ccdt(1,3:end));
        timespectrum(:,1:timespectrum_leng)=datasetfile.dataset.ccdt(place:end,3:end)+timespectrum(:,1:timespectrum_leng);
        timespectrum_normalized(:,1:timespectrum_leng)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range')+timespectrum_normalized(:,1:timespectrum_leng);
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
        timelifetime_leng=length(Lifetime);
        timelifetime(1:timelifetime_leng,len_i)=Lifetime;
        %Intensity change with time
        timeintensity_leng=length(datasetfile.dataset.scatterplot.intensity(1,:));
        timeintensity(len_i,1:timeintensity_leng)=datasetfile.dataset.scatterplot.intensity(1,:);
        %dtime change with time
        SecDtimeintensity_leng=length(SecDtime(:,2));
        SecDtimeintensity(1:SecDtimeintensity_leng,len_i)=SecDtime(:,2);
        
%         %E0001 Ratio change with time
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         timeE0001(len_i,:)=E00sum./E01sum;
    catch ME
        rethrow(ME)
        %disp([name 'may not have length 99'])
    end
end
%%
intensityedge=min(timeintensity(:)):(max(timeintensity(:))-min(timeintensity(:)))/100:max(timeintensity(:));
int10000=1:100:10000;
timeavehis=zeros(max_secs,220);timemaxhis=zeros(max_secs,220);timelifetimehis=zeros(max_secs,(2500-50)/10);timeintensityhis=zeros(max_secs,100);
intSecDtime=zeros(max_secs,6251);
% timeE0001his=zeros(99,20);

%histogramize
for i=1:1:max_secs
    try
        timeavehis(i,:)=histcounts(timeave(i,:),edges);
        timemaxhis(i,:)=histcounts(timemax(:,i),edges);
        timelifetimehis(i,:)=histcounts(timelifetime(i,:),50:10:2500);
        timeintensityhis(i,:)=histcounts(timeintensity(:,i),intensityedge);
        timeint10000his(i,:)=histcounts(timeintensity(:,i),int10000);
        intSecDtime(i,:)=sum(cell2mat(SecDtimeintensity(i,:)'),1);
        % timeE0001his(i,:)=histcounts(timeE0001(:,i),0:0.05:1);
    catch ME
        disp(ME.identifier)
    end
end

meanintensity=mean(timeintensity,1);
stdintensity=std(timeintensity,1);

try
    cd([srdir '/time change/']);
catch
    mkdir([srdir '/time change/']);
    cd([srdir '/time change/']);
end

figure('Position',[680,630,911,348])
subplot(1,2,1)
  surf(edges(1,2:end),1:1:max_secs,timeavehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Average wavelength change with time ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),1:1:max_secs,normalize(timeavehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized average wavelength change with time ' solvent])  
saveas(gcf,[solvent ' average wavelength change with time.jpg']);
  saveas(gcf,[solvent ' average wavelength change with time.fig']);
  disp('Save average wavelength change with time successfully /n');
  close all

figure('Position',[680,630,911,348])
subplot(1,2,1)
  surf(edges(1,2:end),1:1:max_secs,timemaxhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Maximum Wavelength change with time ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),1:1:max_secs,normalize(timemaxhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized maximum wavelength change with time ' solvent])
saveas(gcf,[solvent ' Normalized Maximum change with time.jpg']);
  saveas(gcf,[solvent ' Normalized Maximum change with time.fig']);
  disp('Save Normalized Maximum change with time successfully /n');
  close all
  
figure('Position',[680,630,911,348])
subplot(1,2,1)
  surf(60:10:2500,1:1:max_secs,timelifetimehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Lifeitme change with time ' solvent])
subplot(1,2,2)
  surf(60:10:2500,1:1:max_secs,normalize(timelifetimehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifeitme change with time ' solvent])
saveas(gcf,[solvent ' lifetime change with time.jpg']);
  saveas(gcf,[solvent ' lifetime change with time.fig']);
  disp('Save lifeitme change with time successfully /n');
  close all
  
figure('Position',[680,630,911,348])
subplot(1,2,1)
  surf(intensityedge(1,2:end),1:1:max_secs,timeintensityhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Intensity change with time ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,2:end),1:1:max_secs,normalize(timeintensityhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
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
  
figure('Position',[680,630,911,348])
subplot(1,2,1)
  surf(1:1:max_secs+1,datasetfile.dataset.ccdt(place:end,1),[timespectrum(:,:),zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with time ' solvent])
subplot(1,2,2)
  surf(1:1:max_secs+1,datasetfile.dataset.ccdt(place:end,1),normalize([timespectrum(:,:),zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with time ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with time.fig']);
  disp('Save Spectrum (add up) change with time successfully /n');
  close all  
  
figure('Position',[680,630,911,348])
subplot(1,2,1)
  surf(1:1:max_secs+1,datasetfile.dataset.ccdt(place:end,1),[timespectrum_normalized(:,:),zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with time ' solvent])
subplot(1,2,2)
  surf(1:1:max_secs+1,datasetfile.dataset.ccdt(place:end,1),normalize([timespectrum_normalized(:,:),zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with time ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with time.fig']);
  disp('Save Spectrum (normalize then add up) change with time successfully /n');
  close all   

 
figure
subplot(1,2,1)
  surf(int10000(1,2:end),1:1:max_secs,timeint10000his,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Intensity 10000 change with time ' solvent])
subplot(1,2,2)
  surf(int10000(1,2:end),1:1:max_secs,normalize(timeint10000his,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
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
  
  figure
subplot(1,2,1)
  surf((1:6251)*8/1000,1:max_secs+1,[intSecDtime;zeros(1,6251)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  xlim([0 8])
  title(['Lifetime curve change with int ' solvent])
subplot(1,2,2)
  surf((1:6251)*8/1000,1:max_secs+1,normalize([intSecDtime;zeros(1,6251)],2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  xlim([0 8])
  title(['Normalized lifetime curve change with time ' solvent])
saveas(gcf,[solvent ' Normalized lifetime curve change with time.jpg']);
  saveas(gcf,[solvent ' Normalized lifetime curve change with time.fig']);
  disp('Normalized lifetime curve change with time successfully /n');
  close all   
