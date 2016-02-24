function [y_order, x_order] = forder(XS, Nx)
y_neg = [XS( (Nx/2)+1:Nx)];
y_pos = [XS(1:(Nx/2))];
clear y_order;
y_order = [y_neg, y_pos];
x_neg = -Nx/2+1:1:0;
x_pos = 1:1:Nx/2;
clear x_order;
x_order = [x_neg, x_pos];
return; 


