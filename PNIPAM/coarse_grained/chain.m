function [Cc,CCN] = chain(DP,unimer_C,nC,lambda,a,Cr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% random walk
n=zeros(3,DP);
l=zeros(3,DP);
lc=zeros(3,DP);
P=zeros(3,3);
Y=zeros(3,3);
O=zeros(3,3,DP);
for i=1:DP
    if i==1
        n(:,i) = [1;0;0];
        l(:,i) = n(:,i);
        P(:,:,i) = eye(3);
        Y(:,:,i) = eye(3);
        O(:,:,i) = P(:,:,i)*Y(:,:,i);
    else
        theta_y = 2*(rand-0.5)*pi/a*(2*DP^2+i^2)/DP^2/2;
        theta_z = 2*(rand-0.5)*pi/a*(2*DP^2+i^2)/DP^2/2;
        P(:,:) = [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
        Y(:,:) = [cos(theta_z) -sin(theta_z) 0; sin(theta_z) cos(theta_z) 0; 0 0 1];
        O(:,:,i) = O(:,:,i-1)*P(:,:)*Y(:,:);
        n(:,i) = P(:,:)*Y(:,:)*n(:,i-1);
        l(:,i) = l(:,i-1) + n(:,i);        
        
    end
end
        lc = l*lambda;
% box on
% axis equal
% plot3(l(1,:),l(2,:),l(3,:),  '-o','Color','b','MarkerSize',5,'MarkerFaceColor','#D9FFFF')
% hold on
% plot3(0,0,0, '-o','Color','k','MarkerSize',20,'MarkerFaceColor','#D9FFFF')

%% map unimer
%C
for j=1:DP
        for k=1:nC
            m_unimer_C(:,k,j) = O(:,:,j)*unimer_C(:,k)+lc(:,j)+[Cr;0;0];
        end
end
Cc = reshape(m_unimer_C,[3,DP*nC]);
% Cc = [Cc Cc(:,end)+O(:,:,end)*[2;0;0]];
CCN = Cc(:,end)+O(:,:,end)*[2;0;0];
% %H
% for j=1:DP
%         for k=1:nH
%             m_unimer_H(:,k,j) = O(:,:,j)*unimer_H(:,k)+lc(:,j)+[Cr;0;0];
%         end
% end
% Hc = reshape(m_unimer_H,[3,DP*nH]);
% %N
% for j=1:DP
%         for k=1:nN
%             m_unimer_N(:,k,j) = O(:,:,j)*unimer_N(:,k)+lc(:,j)+[Cr;0;0];
%         end
% end
% Nc = reshape(m_unimer_N,[3,DP*nN]);
% Nc = [Nc Cc(:,end)+O(:,:,end)*[4;0;0]];
% NCN = Nc(:,end)+O(:,:,end)*[4;0;0];
% %O
% for j=1:DP
%         for k=1:nO
%             m_unimer_O(:,k,j) = O(:,:,j)*unimer_O(:,k)+lc(:,j)+[Cr;0;0];
%         end
% end
% Oc = reshape(m_unimer_O,[3,DP*nO]);
end

