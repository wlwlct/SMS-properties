for m=1:38
figure('Position',[2124,1613,1216,326]);
subplot(1,3,1);plot(datasetfile.dataset.ccdt(:,1),TEST_D_pos_sum(:,m));hold on; xline(TEST_EDGE(1,m));
subplot(1,3,2);plot(datasetfile.dataset.ccdt(:,1),TEST_D_neg_sum(:,m));hold on; xline(TEST_EDGE(1,m));
subplot(1,3,3);plot(datasetfile.dataset.ccdt(:,1),TEST_current_sum(:,m))
hold on; plot(datasetfile.dataset.ccdt(:,1),TEST_next_sum(:,m))
hold on; xline(TEST_EDGE(1,m))
legend('Current','Next',num2str(TEST_EDGE(1,m)))
end

p=ceil((1:38)./5);
%figure;
for i=1:8
   figure%subplot(2,4,i);
   plot(TEST_current_sum(:,p==i))
end

figure;surf(normalize(TEST_current_sum,1,'range'),'EdgeColor','none');view([0 0 1]);colormap(jet)
figure;contourf(normalize(TEST_current_sum,1,'range'),[0,0.7,0.9]);colormap(jet)