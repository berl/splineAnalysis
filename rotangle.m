function outmatrix = rotangle(inputangle)

%%%convert angle to radians
inangle = inputangle*2*pi/360;


outmatrix = [ cos(inangle) -sin(inangle); sin(inangle) cos(inangle)];

