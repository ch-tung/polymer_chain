clear;
close all
%% backbone
% Coordinate of C atoms in each unit
% unit_C = load('b_c.dat')';
unit_C = zeros(3,1); % coordinate of C atoms in each unit

% Degree of polymerization
DP_backbone = 100; 

% Chain stiffness
a_backbone = 2;

% Unit persistence
% lambda = 2.5234667395;
lambda = 1;

[lc_backbone,Cc_backbone,O_backbone,n_backbone] = chain(DP_backbone,a_backbone,lambda,unit_C);

% plot backbone
pC_b=plot3(Cc_backbone(1,:),Cc_backbone(2,:),Cc_backbone(3,:), 'o-','Color','#303030','MarkerSize',8,'MarkerFaceColor','#303030','LineWidth',2);

set (gcf,'Position',[0,0,800,800])
axis equal

%% branch
inteval = 1;
DP_branch = 1;
a_branch = 10;

N_branch = floor(DP_backbone/inteval);
Cc_branch = zeros(3,DP_branch,N_branch);
Cc_branch_r = zeros(3,DP_branch,N_branch);

hold on

for i = 1:N_branch
Cc_branch(:,:,i) = chain(DP_branch,a_branch,lambda,unit_C);
Cc_branch(:,:,i) = Cc_branch(:,:,i) - Cc_branch(:,1,i);

t = pi/2;
p = pi*2/7*i;

Y_branch = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
R_branch = [1 0 0; 0 cos(p) sin(p); 0 -sin(p) cos(p)];

Cc_branch_r(:,:,i) = O_backbone(:,:,i*inteval)*R_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
plot3(Cc_branch_r(1,:,i),Cc_branch_r(2,:,i),Cc_branch_r(3,:,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
plot3(Cc_branch_r(1,end,i),Cc_branch_r(2,end,i),Cc_branch_r(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
end

% for i = 1:N_branch
% Cc_branch(:,:,i) = WLM_chain(DP_branch,a_branch,lambda,unit_C);
% Cc_branch(:,:,i) = Cc_branch(:,:,i) - Cc_branch(:,1,i);
% 
% t = pi/12;
% p1 = pi*6/3+pi*2/17*i;
% p2 = pi*5/3+pi*2/17*i;
% p3 = pi*4/3+pi*2/17*i;
% p4 = pi*3/3+pi*2/17*i;
% p5 = pi*2/3+pi*2/17*i;
% p6 = pi*1/3+pi*2/17*i;
% 
% R1_branch = [1 0 0; 0 cos(p1) sin(p1); 0 -sin(p1) cos(p1)];
% R2_branch = [1 0 0; 0 cos(p2) sin(p2); 0 -sin(p2) cos(p2)];
% R3_branch = [1 0 0; 0 cos(p3) sin(p3); 0 -sin(p3) cos(p3)];
% R4_branch = [1 0 0; 0 cos(p4) sin(p4); 0 -sin(p4) cos(p4)];
% R5_branch = [1 0 0; 0 cos(p5) sin(p5); 0 -sin(p5) cos(p5)];
% R6_branch = [1 0 0; 0 cos(p6) sin(p6); 0 -sin(p6) cos(p6)];
% 
% Y_branch = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
% 
% Cc_branch_r1(:,:,i) = O_backbone(:,:,i*inteval)*R1_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
% Cc_branch_r2(:,:,i) = O_backbone(:,:,i*inteval)*R2_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
% Cc_branch_r3(:,:,i) = O_backbone(:,:,i*inteval)*R3_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
% Cc_branch_r4(:,:,i) = O_backbone(:,:,i*inteval)*R4_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
% Cc_branch_r5(:,:,i) = O_backbone(:,:,i*inteval)*R5_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
% Cc_branch_r6(:,:,i) = O_backbone(:,:,i*inteval)*R6_branch*Y_branch*Cc_branch(:,:,i) + lc_backbone(:,i*inteval);
% 
% plot3(Cc_branch_r1(1,2:end,i),Cc_branch_r1(2,2:end,i),Cc_branch_r1(3,2:end,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
% plot3(Cc_branch_r2(1,2:end,i),Cc_branch_r2(2,2:end,i),Cc_branch_r2(3,2:end,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
% plot3(Cc_branch_r3(1,2:end,i),Cc_branch_r3(2,2:end,i),Cc_branch_r3(3,2:end,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
% plot3(Cc_branch_r4(1,2:end,i),Cc_branch_r4(2,2:end,i),Cc_branch_r4(3,2:end,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
% plot3(Cc_branch_r5(1,2:end,i),Cc_branch_r5(2,2:end,i),Cc_branch_r5(3,2:end,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
% plot3(Cc_branch_r6(1,2:end,i),Cc_branch_r6(2,2:end,i),Cc_branch_r6(3,2:end,i), 'o-','Color','#002030','MarkerSize',6,'MarkerFaceColor','#0080C0');
% 
% plot3(Cc_branch_r1(1,end,i),Cc_branch_r1(2,end,i),Cc_branch_r1(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
% plot3(Cc_branch_r2(1,end,i),Cc_branch_r2(2,end,i),Cc_branch_r2(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
% plot3(Cc_branch_r3(1,end,i),Cc_branch_r3(2,end,i),Cc_branch_r3(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
% plot3(Cc_branch_r4(1,end,i),Cc_branch_r4(2,end,i),Cc_branch_r4(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
% plot3(Cc_branch_r5(1,end,i),Cc_branch_r5(2,end,i),Cc_branch_r5(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
% plot3(Cc_branch_r6(1,end,i),Cc_branch_r6(2,end,i),Cc_branch_r6(3,end,i), 'o-','Color','#800000','MarkerSize',8,'MarkerFaceColor','#D00000');
% end

box off
axis off
% xlim([0 100])
% ylim([-60 60])
% zlim([-60 60])
view([0 1 1])
camproj('perspective')
% axis normal;