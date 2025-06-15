function [theta1_1, theta2_1] = IKin2Rfunction(a1, a2, px, py)
c2=(px*px + py*py - a1*a1 - a2*a2)/(2*a1*a2);
%if c2>1, px and py lies outside the workspace of this robot.
if c2<= 1
    s2_1=sqrt(1-c2*c2);
    theta2_1=atan2(s2_1,c2); %first solution
    denom_1=a1*a1 + a2*a2 + 2*a1*a2*cos(theta2_1);

    s1_1= py*(a1+a2*cos(theta2_1)) - px*a2*sin(theta2_1);
    c1_1= px*(a1+a2*cos(theta2_1)) - py*a2*sin(theta2_1);

    theta1_1 = atan2(s1_1, c1_1); %first solution
end
