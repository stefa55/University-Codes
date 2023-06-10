% AEROSPACE DESIGN PROJECT 4 - ANTARCTIC AIRCRAFT
% CENTRE OF GRAVITY AND MOMENTS OF INERTIA CALCULATOR
% Written by S. Messina 2477336M
% James Watt School of Engineering
% University of Glasgow
% - -------------------------------------------------- --
% - -------------------------------------------------- --
% - -------------------------------------------------- --
clear all
clc
close all
% Display user message
fprintf ('CENTRE OF GRAVITY AND MOMENTS OF INERTIA CALCULATOR \n') ;
fprintf (' --------------------------------------------------- \n\n') ;

% M [ mass1 x- pos1 y- pos1 z- pos1 ; mass2 x- pos2 y- pos2 z- pos2 ; etc ...]
M = input ([ 'Please enter all the masses involved using the format \n'...
    '[ mass1 x- pos1 y- pos1 z- pos1 ; mass2 x- pos2 y- pos2 z- pos2 ; etc ...]:\n']) ;
fprintf ('\n') ;
s = size ( M ) ;
mass = M (: ,1);
x_c = M (: ,2) ;
y_c = M (: ,3) ;
z_c = M (: ,4) ;

x_cg = sum ( x_c .* mass ) / sum ( mass ) ;
y_cg = sum ( y_c .* mass ) / sum ( mass ) ;
z_cg = sum ( z_c .* mass ) / sum ( mass ) ;
I_xx = sum ((( y_c - y_cg ) .^2 + ( z_c - z_cg ) .^2) .* mass ) ;
I_yy = sum ((( x_c - x_cg ) .^2 + ( z_c - z_cg ) .^2) .* mass ) ;
I_zz = sum ((( x_c - x_cg ) .^2 + ( y_c - y_cg ) .^2) .* mass ) ;
I_xy = sum (( x_c - x_cg ) .* ( y_c - y_cg ) .* mass ) ;
I_xz = sum (( x_c - x_cg ) .* ( z_c - z_cg ) .* mass ) ;
I_yz = sum (( y_c - y_cg ) .* ( z_c - z_cg ) .* mass ) ;
I = [ I_xx - I_xy - I_xz ; - I_xy I_yy I_yz ; - I_xz - I_yz I_zz ];
[V , D ] = eig ( I ) ;


fprintf('The centre of gravity is: (%.4f, %.4f, %.4f)\n \n', x_cg , y_cg , z_cg ) ;
fprintf('The Matrix of Inertia is :\n \n') ;
disp ( I ) ;
pre = '\it Mass ';
labels = {};
for k = 1: s (1)
labels = [ labels ;[ pre , num2str(k ,'%2d') ]];
end
% labels = [' Aircraft Mass ', 'Right Ski ', 'Left Ski ', 'Cargo Door ', '
% Ramps ', 'Front Ski ', 'Left Door ', 'Right Door '];
% TOP VIEW
img = imread ('BAE146 - Drawing .png ') ;
image ('CData ',img ,'XData ' ,[ -14.39 14.169] , 'YData ' ,[ -13.17 13.17])
hold on
plot ( M (: ,2) , M (: ,3) , 'r.', ' MarkerSize ' ,18) ;
plot ( x_cg , y_cg , 'mp ', ' MarkerSize ', 6 , ' LineWidth ' ,2) ;
xlabel ('\it x - coordinate (m)') ;
ylabel ('\it y - coordinate (m)')
title ('Extra Mass Disposition on the BAe -146 (Top View )', ' FontSize ',14)
% text (M(: ,2) ,M(: ,3) ,labels ,' VerticalAlignment ',' middle ','
% HorizontalAlignment ' ,...
% 'left ', " FontSize " , 10 ," Color " , 'b', 'FontWeight ','bold ');
cog = 'CoG ';
text ( x_cg , y_cg , cog ,' VerticalAlignment ','bottom ',' HorizontalAlignment ',...
    'left ', " FontSize " , 14 ," Color " , 'm', ' FontWeight ','bold ') ;
hold off

% SIDE VIEW
figure
img = imread ('BAE146 - SideM .png ') ;
image ('CData ',img ,'XData ' ,[ -14.39 14.169] , 'YData ' ,[6.31 -2.3])
hold on
plot ( M (: ,2) , M (: ,4) , 'r.', ' MarkerSize ' ,18) ;
xlabel ('\it x - coordinate (m)') ;
ylabel ('\it z - coordinate (m)')
yline (0 , ' -.', ' LineWidth ', 1.2) ;
title ('Extra Mass Disposition on the BAe -146 ( Side View )', ' FontSize ',14)
% text (M(: ,2) ,M(: ,4) ,labels ,' VerticalAlignment ',' middle ','
% HorizontalAlignment ' ,...
% 'left ', " FontSize " , 10 ," Color " , 'b', 'FontWeight ','bold ');
plot ( x_cg , z_cg , 'mp ', ' MarkerSize ', 6 , ' LineWidth ' ,2) ;
cog = 'CoG ';

text ( x_cg , z_cg , cog ,' VerticalAlignment ','bottom ',' HorizontalAlignment ',...
    'left ', " FontSize " , 14 ," Color " , 'm', ' FontWeight ','bold ') ;
hold off

% FRONT VIEW
figure
img = imread ('BAE146 - FrontM .png ') ;
image ('CData ',img ,'XData ' ,[ -13.17 13.17] , 'YData ' ,[6.31 -2.3])
hold on
plot ( M (: ,3) , M (: ,4) , 'r.', ' MarkerSize ' ,18) ;
xlabel ('\it y - coordinate (m)') ;
ylabel ('\it z - coordinate (m)')
yline (0 , ' -.', ' LineWidth ', 1.2) ;
title ('Extra Mass Disposition on the BAe -146 ( Front View )', ' FontSize ',14)
% text (M(: ,3) ,M(: ,4) ,labels ,' VerticalAlignment ',' middle ','
% HorizontalAlignment ' ,...
% 'left ', " FontSize " , 10 ," Color " , 'b', 'FontWeight ','bold ');
plot ( y_cg , z_cg , 'mp ', ' MarkerSize ', 6 , ' LineWidth ', 2) ;
cog = 'CoG ';
text ( y_cg , z_cg , cog ,' VerticalAlignment ','bottom ',' HorizontalAlignment ',...
    'left ', " FontSize " , 14 , " Color " , 'm', ' FontWeight ','bold ');
hold off
