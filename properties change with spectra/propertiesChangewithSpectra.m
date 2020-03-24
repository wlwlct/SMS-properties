%Must have same X% Calculate average wavelength, max wavelength, lifetime
%(choose the first one when calculate in segments), intensity change with
%Intensity. Overall avewav, maxave, lifetime, spectrum(normalized and not normalized.)
%E0001 is not settled

clearvars
solvent='F8T2N2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
srdir=['E:\F8T2400nmCH\'];
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);len=10;

%Find ordered by max; ordered by ratio in three different range.
edges=450:1:670;
year='2020';
place=1;%start to calculate wavelength

Lifindexremove=[];
spectraave=zeros(99,len);
spectramax=zeros(len,99);
spectralifetime=zeros(99,len);
%intE0001=zeros(len,99);
spectraintensity=zeros(len,99);

spectraspectrum=zeros(100-place+1,99,len);
spectraspectrum_normalized=zeros(100-place+1,99,len);
SecDtimeintensity=cell(99,len);
spectraspectrum_diff=zeros(100-place+1,98,len);
spectraspectrum_std=zeros(len,98);
        


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
        %spectraave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
        %max wavelength change with int
        spectramax(len_i,:)=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %spectrum change with int
        spectraspectrum(:,:,len_i)=datasetfile.dataset.ccdt(place:end,3:end);
        spectraspectrum_normalized(:,:,len_i)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range');
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
        spectralifetime(:,len_i)=lifetime;
        %Intensity change with int
        spectraintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
        %dtime change with int
        SecDtimeintensity(:,len_i)=SecDtime(1:99,2);
        %differentce in spectra;jitter
        spectraspectrum_diff(:,:,len_i)=diff(datasetfile.dataset.ccdt(place:end,3:end),1,2);
        spectraspectrum_std(len_i,:)=normalize(std(spectraspectrum_diff(:,:,len_i),0,1),2,'range');
        
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch
        disp([name 'may not have length 99'])
    end
end


spectrum_edge_leng=length(edges)-1;
spectra_max_prepare=cell(spectrum_edge_leng,1);
for spectrum_edge_i=1:spectrum_edge_leng 
    clearvars mol sec
    [mol,sec]=find((spectramax>=edges(1,spectrum_edge_i)) & (spectramax<edges(1,spectrum_edge_i+1)));
    sec_leng=length(sec); spectra_max_prepare{spectrum_edge_i,1}=zeros(100-place+1,sec_leng);
    for sec_i=1:sec_leng
        spectra_max_prepare{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
    end
end

spectra_max_prepare_test=spectra_max_prepare;
for i=1:length(spectra_max_prepare)
    if length(spectra_max_prepare{i,1}(1,:))>=3
        [spectra_max_prepare_test{i,1},~]=sort_spectrum(spectra_max_prepare{i,1},place);
    end
end

T=[];for i=1:220;if ~isempty(spectra_max_prepare{i,1});T=[T,spectra_max_prepare{i,1}];end;end
T_test=[];for i=1:220;if ~isempty(spectra_max_prepare_test{i,1});T_test=[T_test,spectra_max_prepare_test{i,1}];end;end


%edges=430:1:650;
% intE0001his=zeros(99,20);

intensityedge=min(spectraintensity(:)):(max(spectraintensity(:))-min(spectraintensity(:)))/100:max(spectraintensity(:));
intensityedge(1,end)=intensityedge(1,end)+1;%include the last element
intensity_leng=length(intensityedge)-1;
intsp=zeros(100-place+1,intensity_leng);intspn=zeros(100-place+1,intensity_leng);
intavehis=zeros(intensity_leng,220);intmaxhis=zeros(intensity_leng,220);intlifetimehis=zeros(intensity_leng,(2500-50)/10);intintensityhis=zeros(intensity_leng,100);
intSecDtime=zeros(6251,intensity_leng);

for intensity_i=1:intensity_leng
    clearvars mol sec 
    [mol,sec]=find((spectraintensity>=intensityedge(1,intensity_i)) & (spectraintensity<intensityedge(1,intensity_i+1)));
    sec_leng=length(sec);
    intmax_prepare=zeros(sec_leng,1); intave_prepare=zeros(sec_leng,1);intlifetime_prepare=zeros(sec_leng,1);
    for sec_i=1:sec_leng
        intsp(:,intensity_i)=intsp(:,intensity_i)+spectraspectrum(:,sec(sec_i,1),mol(sec_i,1));
        intspn(:,intensity_i)=intspn(:,intensity_i)+spectraspectrum_normalized(:,sec(sec_i,1),mol(sec_i,1));
        intave_prepare(sec_i,1)=spectraave(sec(sec_i,1),mol(sec_i,1));
        intmax_prepare(sec_i,1)=spectramax(mol(sec_i,1),sec(sec_i,1));
        intlifetime_prepare(sec_i,1)=spectralifetime(sec(sec_i,1),mol(sec_i,1));
        intSecDtime(:,intensity_i)=intSecDtime(:,intensity_i)+transpose(SecDtimeintensity{sec(sec_i,1),mol(sec_i,1)});
    end
    intavehis(intensity_i,:)=histcounts(intave_prepare,edges);
    intmaxhis(intensity_i,:)=histcounts(intmax_prepare,edges);
    intlifetimehis(intensity_i,:)=histcounts(intlifetime_prepare,50:10:2500);
end

try
    cd([srdir '/spectra change/']);
catch
    mkdir([srdir '/spectra change/']);
    cd([srdir '/spectra change/']);
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
  histogram(spectraintensity,intensityedge);
  title(['Overall intensity distribution ' solvent])
subplot(2,2,3)
  plot(60:10:2500,sum(intlifetimehis,1),'LineWidth',3);
  title(['Overall lifetime distribution ' solvent]) 
subplot(2,2,4)
  yyaxis left; plot(datasetfile.dataset.ccdt(place:end,1),sum(intsp,2),'LineWidth',3);
  title(['Overall spectrum add up ' solvent])
  yyaxis right;plot(datasetfile.dataset.ccdt(place:end,1),sum(intspn,2),'LineWidth',3);
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
  surf(intensityedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),[intsp(:,:),zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (add up) change with int ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),normalize([intsp(:,:),zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (add up) change with int ' solvent])
saveas(gcf,[solvent ' Spectrum (add up) change with int.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with int.fig']);
  disp('Save Spectrum (add up) change with int successfully /n');
  close all  
  
figure
subplot(1,2,1)
  surf(intensityedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),[intspn(:,:),zeros(100-place+1,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Spectrum (normalize then add up) change with int ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end),datasetfile.dataset.ccdt(place:end,1),normalize([intspn(:,:),zeros(100-place+1,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  title(['Normalized Spectrum (normalize then add up) change with int ' solvent])
saveas(gcf,[solvent ' Spectrum (normalize then add up) change with int.jpg']);
  saveas(gcf,[solvent ' Spectrum (normalize then add up) change with int.fig']);
  disp('Save Spectrum (normalize then add up) change with int successfully /n');
  close all   
 
figure
subplot(1,2,1)
  surf(intensityedge(1,1:end),(1:6251)*8/1000,[intSecDtime zeros(6251,1)],'EdgeColor','none');colormap(jet);view([0 0 1]);
  ylim([0 8])
  title(['Lifetime curve change with int ' solvent])
subplot(1,2,2)
  surf(intensityedge(1,1:end),(1:6251)*8/1000,normalize([intSecDtime zeros(6251,1)],1,'range'),'EdgeColor','none');colormap(jet);view([0 0 1]);
  ylim([0 8])
  title(['Normalized lifetime curve change with int ' solvent])
saveas(gcf,[solvent ' Normalized lifetime curve change with int.jpg']);
  saveas(gcf,[solvent ' Normalized lifetime curve change with int.fig']);
  disp('Normalized lifetime curve change with int successfully /n');
  close all   

  function [Ordered_spectra,Spectra_order]=sort_spectrum(spectraspectrum_normalized,place)
    %compare the distance from each spectrum
    [w,s,m]=size(spectraspectrum_normalized);
    total_spectra=s*m;
    spectra_reshape=reshape(spectraspectrum_normalized,100-place+1,total_spectra);
    spectra_reshape_s=normalize(smoothdata(spectra_reshape,1,'gaussian'),1);
    Spectra_order=[[1;2];zeros(total_spectra-2,1)];
%T=140;figure;plot(spectra_reshape(:,T));hold on;plot(spectra_reshape_test(:,T))
for spectra_order_i=3:total_spectra
    clearvars Point_diff Spectra_order_nonzero Spectra_order_nonzero_leng Point_diff_minloc
    Spectra_order_nonzero=Spectra_order(Spectra_order~=0);
    Spectra_order_nonzero_leng=length(Spectra_order_nonzero);
    Point_diff=sum((spectra_reshape_s(:,spectra_order_i)-spectra_reshape_s(:,Spectra_order_nonzero)).^2,1);
    [~,Point_diff_minloc]=min(Point_diff);
    if Point_diff_minloc==1
        Spectra_order(1:Spectra_order_nonzero_leng+1,1)=[spectra_order_i;Spectra_order_nonzero];
    elseif Point_diff_minloc==Spectra_order_nonzero_leng
        Spectra_order(Spectra_order_nonzero_leng+1,1)=spectra_order_i;
    else
        [~,LorR]=min([Point_diff(1,Point_diff_minloc-1),Point_diff(1,Point_diff_minloc+1)]);
        if LorR==1
            Spectra_order(1:Spectra_order_nonzero_leng+1,1)=[Spectra_order_nonzero(1:Point_diff_minloc-1,1);spectra_order_i;Spectra_order_nonzero(Point_diff_minloc:end,1)];
        elseif LorR==2
            Spectra_order(1:Spectra_order_nonzero_leng+1,1)=[Spectra_order_nonzero(1:Point_diff_minloc,1);spectra_order_i;Spectra_order_nonzero(Point_diff_minloc+1:end,1)];
        else
            disp('Not left or right, something is wrong')
        end
    end     
end
Ordered_spectra=spectra_reshape(:,Spectra_order);
  end