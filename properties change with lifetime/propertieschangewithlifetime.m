%Must have same X% Calculate average wavelength, max wavelength,  intensity change with
%lifetime. %E0001 is not settled

clearvars
solvent='F8T2O2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
srdir=['E:\F8Se2 July\' 'F8Se2O2'];
max_secs=200;%choose according to the max number of seconds in each files.

%srdir=['E:\F8T2O2'];
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);

Lifindexremove=[];
   
lifeave=zeros(max_secs,len);
lifemax=zeros(len,max_secs);
lifelifetime=zeros(max_secs,len);
%intE0001=zeros(len,max_secs);
lifeintensity=zeros(len,max_secs);
place=22;%start to calculate wavelength
lifespectrum=zeros(100-place+1,max_secs,len);
lifespectrum_normalized=zeros(100-place+1,max_secs,len);
edges=450:1:670;

for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    try
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
        lifelifetime_leng=length(lifetime);
        lifelifetime(1:lifelifetime_leng,len_i)=lifetime;
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %average wavelength change with int%each second,each column
        lifeave_leng=length(datasetfile.dataset.scatterplot.spectrum(:,1));
        lifeave(1:lifeave_leng,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with int
        lifemax_leng=length(datasetfile.dataset.ccdt(maxindex+place-1,1));
        lifemax(len_i,1:lifemax_leng)=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %spectrum change with int
        lifespectrum_leng=length(datasetfile.dataset.ccdt(1,3:end));
        lifespectrum(:,1:lifespectrum_leng,len_i)=datasetfile.dataset.ccdt(place:end,3:end);
        lifespectrum_normalized(:,1:lifespectrum_leng,len_i)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range');
        %Intensity change with int
        lifeintensity_leng=length(datasetfile.dataset.scatterplot.intensity(1,:));
        lifeintensity(len_i,1:lifeintensity_leng)=datasetfile.dataset.scatterplot.intensity(1,:);
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch ME
        rethrow(ME)
        disp([name 'may not have length max_secs'])
    end
end

lifetimeedge=min(lifelifetime(:)):(max(lifelifetime(:))-min(lifelifetime(:)))/100:max(lifelifetime(:));
lifetimeedge(1,end)=lifetimeedge(1,end)+1;%include the last element
lifetime_leng=length(lifetimeedge)-1;
lifesp=zeros(100-place+1,lifetime_leng);lifespn=zeros(100-place+1,lifetime_leng);
lifeavehis=zeros(lifetime_leng,220);lifemaxhis=zeros(lifetime_leng,220);
lifeintensityhis=zeros(lifetime_leng,100);
intensity_edge=min(lifeintensity(:)):(max(lifeintensity(:))-min(lifeintensity(:)))/100:max(lifeintensity(:));

for lifetime_i=1:lifetime_leng
    clearvars mol sec 
    [sec,mol]=find((lifelifetime>=lifetimeedge(1,lifetime_i)) & (lifelifetime<lifetimeedge(1,lifetime_i+1)));
    sec_leng=length(sec);
    lifemax_prepare=zeros(sec_leng,1); lifeave_prepare=zeros(sec_leng,1);lifeintensity_prepare=zeros(sec_leng,1);
    for sec_i=1:sec_leng
        lifesp(:,lifetime_i)=lifesp(:,lifetime_i)+lifespectrum(:,sec(sec_i,1),mol(sec_i,1));
        lifespn(:,lifetime_i)=lifespn(:,lifetime_i)+lifespectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
        lifeave_prepare(sec_i,1)=lifeave(sec(sec_i,1),mol(sec_i,1));
        lifemax_prepare(sec_i,1)=lifemax(mol(sec_i,1),sec(sec_i,1));
        lifeintensity_prepare(sec_i,1)=lifeintensity(mol(sec_i,1),sec(sec_i,1));
    end
    lifeavehis(lifetime_i,:)=histcounts(lifeave_prepare,edges);
    lifemaxhis(lifetime_i,:)=histcounts(lifemax_prepare,edges);
    lifeintensityhis(lifetime_i,:)=histcounts(lifeintensity_prepare,intensity_edge);
end

try
    cd([srdir '/lifetime change/']);
catch
    mkdir([srdir '/lifetime change/']);
    cd([srdir '/lifetime change/']);
end

figure
subplot(1,2,1)
  surf(edges(1,2:end),lifetimeedge(1,1:end-1),lifeavehis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['average wavelength change with lifetime ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),lifetimeedge(1,1:end-1),normalize(lifeavehis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized average wavelength change with lifetime ' solvent])  
saveas(gcf,[solvent ' average wavelength change with lifetime.jpg']);
  saveas(gcf,[solvent ' average wavelength change with lifetime.fig']);
  disp('Save average wavelength change with lifetime successfully /n');
  close all

  
figure
subplot(1,2,1)
  surf(edges(1,2:end),lifetimeedge(1,1:end-1),lifemaxhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Maximum Wavelength change with lifetime ' solvent])
subplot(1,2,2)
  surf(edges(1,2:end),lifetimeedge(1,1:end-1),normalize(lifemaxhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Maximum wavelength change with lifetime ' solvent])
saveas(gcf,[solvent ' Normalized Maximum change with lifetime.jpg']);
  saveas(gcf,[solvent ' Normalized Maximum change with lifetime.fig']);
  disp('Save Normalized Maximum change with lifetime successfully /n');
  close all
  
figure
subplot(1,2,1)
  surf(intensity_edge(2:end),lifetimeedge(1,1:end-1),lifeintensityhis,'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['lifeitme change with int ' solvent])
subplot(1,2,2)
  surf(intensity_edge(2:end),lifetimeedge(1,1:end-1),normalize(lifeintensityhis,2,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized lifeitme change with lifetime ' solvent])
saveas(gcf,[solvent ' lifetime change with lifetime.jpg']);
  saveas(gcf,[solvent ' lifetime change with lifetime.fig']);
  disp('Save lifeitme change with lifetime successfully /n');
  close all
    
% figure
% surf(0.05:0.05:1,1:1:max_secs,intE0001his,'EdgeColor','none');
% colormap(jet)
% view([0 0 1])
%   title(['E0001R change with int ' solvent])
%   saveas(gcf,[solvent ' E0001R change with int.jpg']);
%   saveas(gcf,[solvent ' E0001R change with int.fig']);
%   disp('Save E0001R change with int successfully /n');
%   close all  
  
figure
subplot(1,2,1)
  surf(lifetimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),[lifesp(:,:) zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with lifetime ' solvent])
subplot(1,2,2)
  surf(lifetimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),normalize([lifesp(:,:) zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with lifetime ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with lifetime.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with lifetime.fig']);
  disp('Save Spectrum (add up) change with lifetime successfully /n');
  close all  
  
figure
subplot(1,2,1)
  surf(lifetimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),[lifespn(:,:) zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with lifetime ' solvent])
subplot(1,2,2)
  surf(lifetimeedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),normalize([lifespn(:,:) zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with lifetime ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with lifetime.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with lifetime.fig']);
  disp('Save Spectrum (normalize then add up) change with lifetime successfully /n');
  close all   
 
