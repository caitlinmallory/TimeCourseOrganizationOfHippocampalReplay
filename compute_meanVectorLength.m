function [mean_theta,length_r] = compute_meanVectorLength(r,theta)

x = sum(r.*cos(theta))/sum(r);
y = sum(r.*sin(theta))/sum(r);

mean_theta = atan2(y,x);
length_r = sqrt(x^2+y^2);