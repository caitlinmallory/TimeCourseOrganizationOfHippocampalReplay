t_sub = t(t.forward_congruent_replay==1,:);

for i = 1:height(t_sub)
figure()
subplot(2,1,1)
imagesc(t_sub.posterior_left{i}')
caxis([0 0.1])
set(gca,'YDir','normal')
t_sub.ratPos(i)
hold on
yline(2*t_sub.ratPos(i),'r')

hold on
subplot(2,1,2)
imagesc(t_sub.posterior_right{i}')
caxis([0 0.1])
set(gca,'YDir','normal')
t_sub.ratPos(i)
hold on
yline(2*t_sub.ratPos(i),'r')
end


