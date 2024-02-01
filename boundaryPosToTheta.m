function [thetas] = boundaryPosToTheta(x,y,xPos,yPos,dir,a,b)

Y = a.*(-(x - xPos).*sin(dir) + (y - yPos).*cos(dir));
X = b.*((x - xPos).*cos(dir) + (y - yPos).*sin(dir));

thetas = mod(atan2(Y,X),2*pi);