clear
close all

%% Load
data = load('scatter_chain_block_0_g.mat');
data_rod = load('scatter_rod.mat');

q = data_rod.qq;
qL = q*1000;
S_q_rod = data_rod.S_q;
delta_S_q = data_rod.delta_S_q;

S_q_0 = data.S_q;

S_q = zeros(length(q),11*size(S_q_0,2));
p = zeros(3,11*size(S_q_0,2));
for i = 0:10
    filename = ['scatter_chain_block_', num2str(i) ,'_g.mat'];
    data = load(filename);
    S_q(:,1+i*size(S_q_0,2):(i+1)*size(S_q_0,2)) = data.S_q.*delta_S_q';
    p(:,1+i*size(S_q_0,2):(i+1)*size(S_q_0,2)) = data.p;
end

S_q_rod = S_q_rod.*delta_S_q;

%% indexing
b1 = p(2,:)/1000;
rb = p(1,:);
f  = p(3,:);

index_b1 = b1==0.01;
index_rb = rb==10;
index_f = f==0.5;

property = 2;
switch(property)
    case 1
%         index_rb = rb==1;
        index = index_rb&index_f; % b1
        p_color = log10(b1); % b1
    case 2
        index = index_b1&index_f; % rb
        p_color = log10(rb); % rb
    case 3
        index = index_b1&index_rb; % f
        p_color = f; % f
end

i_all = 1:11*size(S_q_0,2);
i_index = i_all(index);

%% color
x_color = (p_color-min(p_color))/(max(p_color)-min(p_color));
color_parula = parula(101);
color = arrayfun(@(x) color_parula(1+round(x*100),:),x_color,'UniformOutput',false);

%% plot S_q
% plot setup
figure(1)
hold on
box on
set(gcf,'Position',[2000,100,600,600])
set(gca,'LineWidth',2)
set(gca,'position',[0.22    0.22   0.72    0.72])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
ylabel('{\it S}({\itQL})','FontSize',24,'Interpreter','tex')
set(gca,'FontSize',28,'FontName','Arial')

% background
for i = 1:9
    q_i = [1e-1 1e3];
    S_q_i_1 = 10^(i-5)./q_i;
    S_q_i_2 = 10^(2*i-10)./q_i.^2;
    plot(q_i,S_q_i_1,'--','LineWidth',1,'Color','#C0C0C0')
    plot(q_i,S_q_i_2,':','LineWidth',2,'Color','#C0C0C0')
end

% S_q_rod
plot(qL,S_q_rod,'k-','LineWidth',2)

% S_q
for i = 1:length(i_index)
    plot(qL,S_q(:,i_index(i)),'-','LineWidth',2,'Color',color{i_index(i)})
end

xlim([1e-1,1e3])
ylim([1e-3,2e0])
xticks([1e-1,1e0,1e1,1e2,1e3])
yticks([1e-3,1e-2,1e-1,1e0])

%% plot delta S_q
% plot setup
figure(2)
hold on
box on
set(gcf, 'color', 'none')
set(gcf,'Position',[2000,100,400,400])
set(gca,'LineWidth',2)
set(gca,'position',[0.25    0.25   0.6    0.6])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
ylabel('\Delta{\itS}({\itQL})','FontSize',24,'Interpreter','tex')
set(gca,'FontSize',24,'FontName','Arial')

% S_q_rod
plot(qL,S_q_rod./S_q_rod,'k-','LineWidth',2)

% S_q
for i = 1:length(i_index)
    plot(qL,S_q(:,i_index(i))./S_q_rod','-','LineWidth',2,'Color',color{i_index(i)})
end

xlim([1e-1,1e3])
ylim([1e0,5e0])
xticks([1e-1,1e0,1e1,1e2,1e3])
% yticks([1e0,1e1])