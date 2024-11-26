function head_direction = calculate_headDirection(left_led_x,left_led_y,right_led_x,right_led_y)

% Calculate head direction

% subtract the position of the left LED so that you have a vector pointing
% from (0,0) in the angle of the red LED.
% compute the angle between a straight line pointing horizontal rightward,
% and the vector going from the left to the right led.

x1 = 1;
y1 = 0;

x2 = right_led_x - left_led_x;
y2 = right_led_y - left_led_y;

angle0 = atan2(x1.*y2-y1.*2,x1.*x2+y1.*y2);
angle0(angle0<0) = angle0(angle0<0)+2*pi; % have all head directions be positive;
% add 90 degrees to get the direction that the rat is heading relative to
% horizontal
head_direction = angle0 + pi/2;
head_direction(head_direction>2*pi) = head_direction(head_direction>2*pi)-2*pi;
