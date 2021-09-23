function [Ct,Ht] = chain_tail(unimer_Ct,unimer_Ht,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% random walk
Ct(1,:) = unimer_Ct(1,:);
Ct(2,:) = unimer_Ct(2,:)+b*cos(unimer_Ct(1,:)/2);
Ct(3,:) = unimer_Ct(2,:)+b*sin(unimer_Ct(1,:)/2);

Ht(1,:) = unimer_Ht(1,:);
Ht(2,:) = unimer_Ht(2,:)+b*cos(unimer_Ht(1,:)/2);
Ht(3,:) = unimer_Ht(2,:)+b*sin(unimer_Ht(1,:)/2);