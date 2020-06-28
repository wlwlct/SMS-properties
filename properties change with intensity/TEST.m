%Must have same X% Calculate average wavelength, max wavelength, lifetime
%(choose the first one when calculate in segments), intensity change with
%Intensity. Overall avewav, maxave, lifetime, spectrum(normalized and not normalized.)
%E0001 is not settled

clearvars
solvent='F8T2N2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\F8T2N2High'];
cd (srdir)
intnum=100;%Determine how many spectra you want to output
allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);

Lifindexremove=[];
   
intave=zeros(99,len);
intmax=zeros(len,99);
intlifetime=zeros(99,len);
%intE0001=zeros(len,99);
intintensity=zeros(len,99);
place=22;%start to calculate wavelength
intspectrum=zeros(100-place+1,99,len);
intspectrum_normalized=zeros(100-place+1,99,len);
SecDtimeintensity=cell(99,len);
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
        [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %average wavelength change with int%each second,each column
        intave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with int
        intmax(len_i,:)=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %spectrum change with int
        intspectrum(:,:,len_i)=datasetfile.dataset.ccdt(place:end,3:end);
        intspectrum_normalized(:,:,len_i)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range');
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
        %intintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
        intintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(2,:);
        %dtime change with int
        SecDtimeintensity(:,len_i)=SecDtime(1:99,2);
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch
        disp([name 'may not have length 99'])
    end
end

save([solvent ' correlation.mat'],'intspectrum','intspectrum_normalized','intintensity','SecDtimeintensity');
%edges=430:1:650;
% intE0001his=zeros(99,20);

intsint=sort(intintensity(:));
clearvars mol sec
[mol,sec]=find(intintensity>=intsint(end-intnum,1));
sec_leng=length(sec);
intsp=zeros(100-place+1,sec_leng);intspn=zeros(100-place+1,sec_leng);
intSecDtime=zeros(6251,sec_leng);intsint=zeros(1,sec_leng);
for sec_i=1:sec_leng
    intsint(1,sec_i)=intintensity(mol(sec_i,1),sec(sec_i,1));
    intsp(:,sec_i)=intspectrum(:,sec(sec_i,1),mol(sec_i,1));
    intspn(:,sec_i)=intspectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
    intSecDtime(:,sec_i)=transpose(SecDtimeintensity{sec(sec_i,1),mol(sec_i,1)});
end
[intsint,intasidx]=sort(intsint);
intsp=intsp(:,intasidx);
intspn=intspn(:,intasidx);
intSecDtime=intSecDtime(:,intasidx);

try
    cd([srdir '/intensity change/']);
catch
    mkdir([srdir '/intensity change/']);
    cd([srdir '/intensity change/']);
end
  
  
figure
subplot(1,2,1)
  surf(1:(sec_leng+1),datasetfile.dataset.ccdt(place:end,1),[intsp(:,:),zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) with large int ' solvent])
subplot(1,2,2)
  surf(1:(sec_leng+1),datasetfile.dataset.ccdt(place:end,1),normalize([intsp(:,:),zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) with large int ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) with large int.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) with large int.fig']);
  disp('Save Spectrum (add up) with large int successfully /n');
  close all  
   
figure
subplot(1,2,1)
  surf(1:(sec_leng+1),(1:6251)*8/1000,[intSecDtime zeros(6251,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  ylim([0 8])
  title(['Lifetime curve with large int ' solvent])
subplot(1,2,2)
  surf(1:(sec_leng+1),(1:6251)*8/1000,normalize([intSecDtime zeros(6251,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  ylim([0 8])
  title(['Normalized lifetime curve change with int ' solvent])
saveas(gcf,[solvent ' Normalized lifetime curve with large int.jpg']);
  saveas(gcf,[solvent ' Normalized lifetime curve with large int.fig']);
  disp('Normalized lifetime curve with large int successfully /n');
  close all   
