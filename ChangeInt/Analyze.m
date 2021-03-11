clearvars -except molecules_CND wl
codefolder='C:\Users\Livi\Documents\GitHub\Data-Reform\After Dataset compare Excel';
edges=450:1:670;
Folder='E:\MEH substrate clean mat data\MEH_Chloroform_rmBG\change int';
% Plot_inc_dec(molecules_CND,wl,edges,codefolder)

wl_leng=length(wl(:,1));
CND_leng=length(molecules_CND(:,1));  total_leng=0;
for CND_leng_i=1:CND_leng;total_leng=total_leng+length(molecules_CND{CND_leng_i,2}(1,:));end

molecules_Diff=zeros(wl_leng,total_leng);
molecules_current=zeros(wl_leng,total_leng);
molecules_next=zeros(wl_leng,total_leng);

current_column=1;
for CND_leng_i=1:CND_leng
    clearvars current_leng
    current_leng=length(molecules_CND{CND_leng_i,2}(1,:));
    molecules_Diff(:,current_column:current_column+current_leng-1)=molecules_CND{CND_leng_i,2};
    molecules_current(:,current_column:current_column+current_leng-1)=molecules_CND{CND_leng_i,1}(:,1:current_leng);
    molecules_next(:,current_column:current_column+current_leng-1)=molecules_CND{CND_leng_i,1}(:,2:current_leng+1);
    current_column=current_column+current_leng;
end

molecules_Diff=rmLowInt(molecules_Diff,codefolder);
s_diff=sum(molecules_Diff,1);
% s_diff_sort=zeros(2,181);
% [s_diff_sort(1,:),s_diff_sort(2,:)]=sort(s_diff);
% figure;mesh(normalize(molecules_Diff(:,s_diff_sort(2,:)),1,'range'));view([0 0 1]);colormap(jet);
decrease_loc=find(s_diff<0);
increase_loc=find(s_diff>=0);

molecules_Diff_decrease=molecules_Diff(:,decrease_loc)*(-1);
molecules_current_decrease=molecules_current(:,decrease_loc);
molecules_next_decrease=molecules_next(:,decrease_loc);

molecules_Diff_increase=molecules_Diff(:,increase_loc);
molecules_current_increase=molecules_current(:,increase_loc);
molecules_next_increase=molecules_next(:,increase_loc);

%Increase calculation order by diff
[increase_diff_prepare,increase_current_prepare,increase_next_prepare]=Spectra_prepare(molecules_Diff_increase,molecules_current_increase,molecules_next_increase,wl,edges);
[increase_diff_mesh,increase_xl]=PreparePlot(increase_diff_prepare,edges,wl);
increase_current_mesh=PreparePlot(increase_current_prepare,edges,wl);
increase_next_mesh=PreparePlot(increase_next_prepare,edges,wl);
increase_xl=cellfun(@num2str,num2cell(increase_xl),'UniformOutput',false);
%Decrease calculation order by diff
[decrease_diff_prepare,decrease_current_prepare,decrease_next_prepare]=Spectra_prepare(molecules_Diff_decrease,molecules_current_decrease,molecules_next_decrease,wl,edges);
[decrease_diff_mesh,decrease_xl]=PreparePlot(decrease_diff_prepare,edges,wl);
decrease_current_mesh=PreparePlot(decrease_current_prepare,edges,wl);
decrease_next_mesh=PreparePlot(decrease_next_prepare,edges,wl);
decrease_xl=cellfun(@num2str,num2cell(decrease_xl),'UniformOutput',false);


%figure('Position',[2582,1003,1440,385]);
figure;
subplot(1,3,1);mesh(1:length(decrease_current_mesh(1,:)),wl,normalize(decrease_current_mesh,1,'range'));view([0 0 1]);colormap(jet);title('Decrease Current');
ylabel('Wavelength (nm)');xticks(1:10:length(decrease_current_mesh(1,:)));xticklabels(decrease_xl(1:10:end));
subplot(1,3,2);mesh(1:length(decrease_current_mesh(1,:)),wl,normalize(decrease_next_mesh,1,'range'));view([0 0 1]);colormap(jet);title('Decrease Next');
ylabel('Wavelength (nm)');xticks(1:10:length(decrease_current_mesh(1,:)));xticklabels(decrease_xl(1:10:end));
subplot(1,3,3);mesh(1:length(decrease_current_mesh(1,:)),wl,normalize(decrease_diff_mesh,1,'range'));view([0 0 1]);colormap(jet);title('Decrease Diff');
ylabel('Wavelength (nm)');xticks(1:10:length(decrease_current_mesh(1,:)));xticklabels(decrease_xl(1:10:end));
% title('decrease')

%figure('Position',[2582,1003,1440,385]);
figure;
subplot(1,3,1);mesh(1:length(increase_current_mesh(1,:)),wl,normalize(increase_current_mesh,1,'range'));view([0 0 1]);colormap(jet);title('increase Current');
ylabel('Wavelength (nm)');xticks(1:10:length(increase_current_mesh(1,:)));xticklabels(increase_xl(1:10:end));
subplot(1,3,2);mesh(1:length(increase_current_mesh(1,:)),wl,normalize(increase_next_mesh,1,'range'));view([0 0 1]);colormap(jet);title('increase Next');
ylabel('Wavelength (nm)');xticks(1:10:length(increase_current_mesh(1,:)));xticklabels(increase_xl(1:10:end));
subplot(1,3,3);mesh(1:length(increase_current_mesh(1,:)),wl,normalize(increase_diff_mesh,1,'range'));view([0 0 1]);colormap(jet);title('increase Diff');
ylabel('Wavelength (nm)');xticks(1:10:length(increase_current_mesh(1,:)));xticklabels(increase_xl(1:10:end));
% title('increase')

% %Increase calculation order by diff
% molecules_Diff_increase=rmLowInt(molecules_Diff_increase,codefolder);
% [increase_diff_prepare,increase_current_prepare,increase_next_prepare]=Spectra_prepare(molecules_Diff_increase,molecules_current_increase,molecules_next_increase,wl,edges);
% increase_diff_mesh=PreparePlot(increase_diff_prepare,edges,wl);
% increase_current_mesh=PreparePlot(increase_current_prepare,edges,wl);
% increase_next_mesh=PreparePlot(increase_next_prepare,edges,wl);
% %Decrease calculation order by diff
% molecules_Diff_decrease=rmLowInt(molecules_Diff_decrease,codefolder);
% [decrease_diff_prepare,decrease_current_prepare,decrease_next_prepare]=Spectra_prepare(molecules_Diff_decrease,molecules_current_decrease,molecules_next_decrease,wl,edges);
% decrease_diff_mesh=PreparePlot(decrease_diff_prepare,edges,wl);
% decrease_current_mesh=PreparePlot(decrease_current_prepare,edges,wl);
% decrease_next_mesh=PreparePlot(decrease_next_prepare,edges,wl);

%%
%plot everything in the decrease and increase mesh
figure;
sta='decrease';
eval(['M=' sta '_diff_mesh;']);
eval(['l=' sta '_xl;']);
mesh_leng=length(M(1,:));
% for i=1:1:mesh_leng;
%     hold on;plot(wl,normalize(M(:,i),1,'range'),'LineWidth',3,'DisplayName',l{1,i});
% end;

close all
flag=1;
for mesh_leng_i=1:mesh_leng
    ceil_num=ceil(mesh_leng_i/4);
    if flag~=ceil_num
        flag=ceil_num;
        legend;xlabel('Wavelength (nm)');ylabel('Normalized Intensity')
        cd(Folder)
        %saveas(gcf,['MEH CH Clear' ' current spectra' num2str(flag) '.jpg']);
        %saveas(gcf,['MEH CH Clear' ' current spectra' num2str(flag) '.fig']);
        close gcf
    end
    figure(ceil_num)
    hold on;plot(wl,normalize(M(:,mesh_leng_i),1,'range'),'LineWidth',3,'DisplayName',l{1,mesh_leng_i});
    title(['MEH CH Clear ' sta 'spectra sort by wavelength at peak maxima'])
end


%%

function FE(molecules_CND)
CND_leng=length(molecules_CND(:,1)); 
First_end=cell(CND_leng,1);
for CND_leng_i=1:CND_leng
    if length(molecules_CND{CND_leng_i,1}(1,:))>1
        First_end{CND_leng_i,1}=molecules_CND{CND_leng_i,1}(:,[1 end]);
    else
        First_end{CND_leng_i,1}=molecules_CND{CND_leng_i,1}(:,1);
    end
    %figure;plot(normalize(First_end{CND_leng_i,1}(:,1)));hold on;plot(normalize(First_end{CND_leng_i,1}(:,2)));
end
end


function spectra1=rmLowInt(spectra,codefolder)
    [spectra_zong,spectra_heng]=size(spectra);
    spectra_stage=zeros(spectra_zong,spectra_heng);
    spectra_stage_ratio=zeros(1,spectra_heng);
    for i=1:spectra_heng
        cd(codefolder)
        clearvars A
        [~,A.eff_fit,~,A.numst,~]=Traceson(spectra(:,i),codefolder);
        if ~isempty(A)
            if length(A.eff_fit(:,1))<A.numst;efffit=A.eff_fit(1,:);else;efffit=A.eff_fit(A.numst,:);end
            spectra_stage(:,i)=transpose(efffit);
            spectra_stage_ratio(1,i)=max(spectra_stage(:,i)+abs(min(spectra_stage(:,i))));
        end
    end
    
    max_int=max(spectra(:));
    spectra(end-5:end,spectra_stage_ratio<1.3)=1000+max_int;
    spectra1=spectra;
end
function [spc_diff_prepare,spc_current_prepare,spc_next_prepare]=Spectra_prepare(Spectra_Diff,Spectra_current,Spectra_Next,wl,edges)
%change of the spectrum, increase or decrease; plot current and next along
%with difference of spectrum
    [Spectra_Diff_zong,~]=size(Spectra_Diff);
    [~,spectramax_smooth_loc]=max(smoothdata(Spectra_Diff,1,'gaussian',8),[],1);
    spc_max_smooth=transpose(wl(spectramax_smooth_loc,1));
%Order by the change of the diff of the spectrum
    spectra_edge_leng=length(edges)-1;
    spc_diff_prepare=cell(spectra_edge_leng,1);
    spc_current_prepare=cell(spectra_edge_leng,1);
    spc_next_prepare=cell(spectra_edge_leng,1);
    for spectra_edge_i=1:spectra_edge_leng 
        clearvars spc
        spc=find((spc_max_smooth>=edges(1,spectra_edge_i)) & (spc_max_smooth<edges(1,spectra_edge_i+1)));
        spc_leng=length(spc); 
        spc_diff_prepare{spectra_edge_i,1}=zeros(Spectra_Diff_zong,spc_leng);
        spc_current_prepare{spectra_edge_i,1}=zeros(Spectra_Diff_zong,spc_leng);
        spc_next_prepare{spectra_edge_i,1}=zeros(Spectra_Diff_zong,spc_leng);
        for sec_i=1:spc_leng
            spc_diff_prepare{spectra_edge_i,1}=Spectra_Diff(:,spc(1,sec_i));
            spc_current_prepare{spectra_edge_i,1}=Spectra_current(:,spc(1,sec_i));
            spc_next_prepare{spectra_edge_i,1}=Spectra_Next(:,spc(1,sec_i));
        end
    end
end
function [mesh_prepare,labs]=PreparePlot(prepare_cell,edges,wl)
    edges_leng=length(edges)-1;
    total_spc=0;
    for edge_leng_i=1:edges_leng
        if ~isempty(prepare_cell{edge_leng_i,1})
        total_spc=total_spc+length(prepare_cell{edge_leng_i,1}(1,:));
        end
    end
    mesh_prepare=zeros(length(wl),total_spc);
    labs=zeros(1,total_spc);
    current_i=1;
    for edge_leng_i=1:edges_leng
        if ~isempty(prepare_cell{edge_leng_i,1})
            current_leng=length(prepare_cell{edge_leng_i,1}(1,:));
            labs(1,current_i:current_i+current_leng-1)=edges(1,edge_leng_i);
            mesh_prepare(:,current_i:current_i+current_leng-1)=prepare_cell{edge_leng_i,1};
            current_i=current_i+current_leng;
        end
    end
end
