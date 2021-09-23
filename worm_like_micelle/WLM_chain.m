function [lc,Cc,O,n] = WLM_chain(DP,a,lambda,unit_C)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
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
        R(:,:,i) = eye(3);
        O(:,:,i) = P(:,:,i)*Y(:,:,i);
    else
        theta_y = 2*(rand-0.5)*pi/a/2;
        theta_z = 2*(rand-0.5)*pi/a/2;
        phi = 2*(rand-0.5)*pi/5;

        P(:,:) = [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
        Y(:,:) = [cos(theta_z) -sin(theta_z) 0; sin(theta_z) cos(theta_z) 0; 0 0 1];
        R(:,:) = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
        
        O(:,:,i) = P*Y*O(:,:,i-1);
        n(:,i) = P*Y*n(:,i-1);
        l(:,i) = l(:,i-1) + n(:,i);
        
    end
end
lc = l*lambda;
%% map unimer
%C
nC = size(unit_C,2);
for j=1:DP
        for k=1:nC
            m_backbone_C(:,k,j) = O(:,:,j)*unit_C(:,k)+lc(:,j)+[0;0;0];
        end
end
Cc = reshape(m_backbone_C,[3,DP*nC]);
end