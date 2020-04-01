%Parameter setting
clearvars;solvent='F8T2N2';
codefolder=pwd;
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
spectramax=zeros(len,99);
spectralifetime=zeros(99,len);
%intE0001=zeros(len,99);
spectraintensity=zeros(len,99);

spectraspectrum=zeros(100-place+1,99,len);
spectraspectrum_normalized=zeros(100-place+1,99,len);
SecDtimeintensity=cell(99,len);
spectraspectrum_diff=zeros(100-place+1,98,len);
spectraspectrum_std=zeros(len,98);
spectra_stage=zeros(100-place+1,99,len);
spectra_stage_ratio=zeros(99,len);

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
        [~,maxindex]=max(datasetfile.dataset.ccdt(place:end,3:end),[],1);
        %max wavelength change with spectra
        spectramax(len_i,:)=datasetfile.dataset.ccdt(maxindex+place-1,1);
        %spectrum change with spectra
        spectraspectrum(:,:,len_i)=datasetfile.dataset.ccdt(place:end,3:end);
        for ii=1:99
            clearvars A
            cd(codefolder)
            [~,A.eff_fit,~,A.numst,~]=Traceson(spectraspectrum(:,ii,len_i),codefolder);
            if ~isempty(A)
                if length(A.eff_fit(:,1))<A.numst;efffit=A.eff_fit(1,:);else;efffit=A.eff_fit(A.numst,:);end
                spectra_stage(:,ii,len_i)=transpose(efffit);
                spectra_stage_ratio(ii,len_i)=max(spectra_stage(:,ii,len_i))/min(spectra_stage(:,ii,len_i));
            end
        end
        spectraspectrum_normalized(:,:,len_i)=normalize(datasetfile.dataset.ccdt(place:end,3:end),1,'range');
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
        spectralifetime(:,len_i)=lifetime;
        %Intensity change with spectra
        spectraintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
        %dtime change with spectra
        %SecDtimeintensity(:,len_i)=SecDtime(1:99,2);
        %differentce in spectra;jitter
        spectraspectrum_diff(:,:,len_i)=diff(datasetfile.dataset.ccdt(place:end,3:end),1,2);
        spectraspectrum_std(len_i,:)=std(spectraspectrum_diff(:,:,len_i),0,1);
%         %E0001 Ratio change with int
%         E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
%         E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
%         intE0001(len_i,:)=E00sum./E01sum;
    catch
        disp([name 'may not have length 99'])
    end
end

%Remove some of the spectrum that has no shape
% spectraspectrum% spectra_stage
% spectraspectrum_diff% spectraspectrum_std
% spectra_stage_ratio% spectraspectrum_normalized% spectramax
% spectralifetime% spectraintensity% SecDtimeintensity

max_int=max(spectraspectrum(:));
[Bad_sec,Bad_mol]=find(spectra_stage_ratio<2.3);
Bad_leng=length(Bad_sec);
for i=1:Bad_leng
    spectraspectrum(end-10:end,Bad_sec(i,1),Bad_mol(i,1))=100+max_int;
    spectraspectrum(end-10:end,Bad_sec(i,1),Bad_mol(i,1))=100+max_int;
    spectraspectrum_normalized(90:100,Bad_sec(i,1),Bad_mol(i,1))=1.1;
end


%cut anything related to the last second, it would not related to spectrum
%change or jitter.
spectramax=spectramax(:,1:98);
spectramax_smooth=spectramax_smooth(:,1:98);
spectraspectrum_normalized_smooth=smoothdata(spectraspectrum_normalized,1,'gaussian',8);
%spectraspectrum=spectraspectrum(:,1:98,:);
%spectraspectrum_normalized=spectraspectrum_normalized(:,1:98,:);
spectralifetime=spectralifetime(1:98,:);
%spectraintensity=spectraintensity(:,1:98);
%SecDtimeintensity=SecDtimeintensity(1:98,:);
spectraspectrum_std_normalize=normalize(spectraspectrum_std,2,'range');

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
        spectra_next_intensity_prepare{spectrum_edge_i}(1,sec_i)=spectraintensity(mol(sec_i,1),sec(sec_i,1)+1);
        spectra_jitter_prepare{spectrum_edge_i}(1,sec_i)=spectraspectrum_std(mol(sec_i,1),sec(sec_i,1));
        spectra_jitter_prepare_normalize{spectrum_edge_i}(1,sec_i)=spectraspectrum_std_normalize(mol(sec_i,1),sec(sec_i,1));
        spectra_next_prepare{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized(:,sec(sec_i,1)+1,mol(sec_i,1));
        spectra_next_prepare_smooth{spectrum_edge_i}(:,sec_i)=spectraspectrum_normalized_smooth(:,sec(sec_i,1)+1,mol(sec_i,1));
        spectra_difference_prepare{spectrum_edge_i}(:,sec_i)=spectraspectrum_diff(:,sec(sec_i,1),mol(sec_i,1));
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
spectra_difference_average=zeros(100-place+1,220);
for i=1:220
    if ~isempty(spectra_max_prepare{i,1})
        spectra_max_average(:,i)=mean(spectra_max_prepare{i,1},2);
        spectra_next_average(:,i)=mean(spectra_next_prepare{i,1},2);
        spectra_difference_average(:,i)=mean(spectra_difference_prepare{i,1},2);
    end
end
figure('Position',[0,0,762,211]);
subplot(1,3,1);mesh(1:220,datasetfile.dataset.ccdt(place:end,1),spectra_max_average);
    view([0 0 1]); colormap(jet);title('current')
subplot(1,3,2);mesh(1:220,datasetfile.dataset.ccdt(place:end,1),spectra_next_average);
    view([0 0 1]); colormap(jet);title('next')
subplot(1,3,3);mesh(1:220,datasetfile.dataset.ccdt(place:end,1),spectra_difference_average);
    view([0 0 1]); colormap(jet);title('difference')
saveas(gcf,['current next and difference spectra.fig']);
saveas(gcf,['current next and difference spectra.jpg']);
close all

%In each max wavelength range, positive, negative, current and next
spectra_D_pos=spectra_difference_prepare;
spectra_D_neg=spectra_difference_prepare;
spectra_D_pos_sum=zeros(100-place+1,220);
spectra_D_neg_sum=zeros(100-place+1,220);
spectra_max_sum=zeros(100-place+1,220);
spectra_next_sum=zeros(100-place+1,220);

for i=1:220
    if ~isempty(spectra_max_prepare{i,1})
        spectra_D_pos{i,1}(spectra_D_pos{i,1}<0)=0;
        spectra_D_neg{i,1}(spectra_D_neg{i,1}>0)=0;
        spectra_D_pos_sum(:,i)=sum(spectra_D_pos{i,1},2);
        spectra_D_neg_sum(:,i)=sum(spectra_D_neg{i,1},2);
        spectra_max_sum(:,i)=sum(spectra_max_prepare{i,1},2);
        spectra_next_sum(:,i)=sum(spectra_next_prepare{i,1},2);
    end
end

figure('Position',[0,0,1032,741]);
subplot(2,2,1);mesh(1:220,datasetfile.dataset.ccdt(:,1),normalize(spectra_D_pos_sum,1,'range'));
view([0 0 1]);colormap(jet);title(['Spectra increase'])
subplot(2,2,2);mesh(1:220,datasetfile.dataset.ccdt(:,1),normalize(spectra_D_neg_sum*(-1),1,'range'));
view([0 0 1]);colormap(jet);title('Spectra decrease')
subplot(2,2,3);mesh(1:220,datasetfile.dataset.ccdt(:,1),normalize(spectra_max_sum,1,'range'));
view([0 0 1]);colormap(jet);title('Current Sum')
subplot(2,2,4);mesh(1:220,datasetfile.dataset.ccdt(:,1),normalize(spectra_next_sum,1,'range'));
view([0 0 1]);colormap(jet);title('Next Sum')

saveas(gcf,['Pos Neg current next.fig']);
saveas(gcf,['Pos Neg current next.jpg']);
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
%%
% Use to sort spectrum;spectra_max_prepare_test=spectra_max_prepare;
% for i=1:length(spectra_max_prepare)
%     if length(spectra_max_prepare{i,1}(1,:))>=3
%         [spectra_max_prepare_test{i,1},~]=sort_spectrum(spectra_max_prepare{i,1},place);
%     end
% end
%%
% T=[];for i=1:220;if ~isempty(spectra_max_prepare{i,1});T=[T,spectra_max_prepare{i,1}];end;end
% T_smooth=[];for i=1:220;if ~isempty(spectra_max_prepare{i,1});T_smooth=[T_smooth,spectra_max_prepare_smooth{i,1}];end;end
% I_NEXT=[];for i=1:220;if ~isempty(spectra_jitter_prepare{i,1});I_NEXT=[I_NEXT,spectra_next_intensity_prepare{i,1}];end;end
% J=[];for i=1:220;if ~isempty(spectra_jitter_prepare{i,1});J=[J,spectra_jitter_prepare{i,1}];end;end
% J_N=[];for i=1:220;if ~isempty(spectra_jitter_prepare_normalize{i,1});J_N=[J_N,spectra_jitter_prepare_normalize{i,1}];end;end
% % T_test=[];for i=1:220;if ~isempty(spectra_max_prepare_test{i,1});T_test=[T_test,spectra_max_prepare_test{i,1}];end;end
% N=[];for i=1:220;if ~isempty(spectra_next_prepare{i,1});N=[N,spectra_next_prepare{i,1}];end;end
% N_smooth=[];for i=1:220;if ~isempty(spectra_next_prepare{i,1});N_smooth=[N_smooth,spectra_next_prepare_smooth{i,1}];end;end
% D=[];for i=1:220;if ~isempty(spectra_difference_prepare{i,1});D=[D,spectra_difference_prepare{i,1}];end;end
% %In max wavelenght range, intensity vs: max wavelength range
% range=5:150;
% TEST_N_smooth=N_smooth(:,range);TEST_T_smooth=T_smooth(:,range);
% TEST_N=N(:,range);TEST_T=T(:,range);TEST_D=D(:,range);
% [~,D_max_idx]=max(TEST_N_smooth,[],1);D_max=datasetfile.dataset.ccdt(D_max_idx+place-1,1);
% D_intensity=sum(TEST_D,1);
% figure;scatter(D_intensity,D_max)
% %plot original spectrum and next spectrum for those max wavelength change
% D_changeloc=find(D_max>520);
% I_changeloc=find(D_intensity>-800);
% changeloc=intersect(D_changeloc,I_changeloc);
% for i=1:length(D_changeloc)
%     figure('Position',[2293,1127,1762,401]);
%     subplot(1,3,1);plot(datasetfile.dataset.ccdt(:,1),TEST_T(:,changeloc(i,1)),'DisplayName',['origin' num2str(D_changeloc(i,1))]);
%     hold on;plot(datasetfile.dataset.ccdt(:,1),TEST_N(:,changeloc(i,1)),'DisplayName',['next' num2str(D_changeloc(i,1))]);
%     legend
%     subplot(1,3,2);plot(datasetfile.dataset.ccdt(:,1),TEST_T_smooth(:,changeloc(i,1)),'DisplayName',['origin' num2str(D_changeloc(i,1))]);
%     hold on;plot(datasetfile.dataset.ccdt(:,1),TEST_N_smooth(:,changeloc(i,1)),'DisplayName',['next' num2str(D_changeloc(i,1))]);
%     legend
%     subplot(1,3,3);plot(datasetfile.dataset.ccdt(:,1),TEST_D(:,changeloc(i,1)));
% end   
%%  
%Directly calculate std of each spectrum
%remove some with very small std% %Give higher weight for 480 and 590
% [~,loc480]=min(abs(datasetfile.dataset.ccdt(:,1)-480));
% [~,loc590]=min(abs(datasetfile.dataset.ccdt(:,1)-590));
% w=ones(length(datasetfile.dataset.ccdt(:,1)),1);
% w(loc480:loc590)=2;
% spectra_std_w=zeros(94,99);
% for i=1:94
%    spectra_std_w(i,:)=std(spectraspectrum(:,:,i),w,1);
% end