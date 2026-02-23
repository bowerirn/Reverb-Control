% This function computes the sound pressure from a panel of size (Lx,Ly) at
% a point in space (x,y,z), where z is the distance outward from the panel
% center.

% F is the frequency at which the calculation takes place (Hz)
% A is the panel displacement matrix 
% (for simplicity, use linspace so panel is 100x100 matrix)

function out = discrete_Rayleigh(A,Lx,Ly,f,x,y,z)

% convert to radians
w = 2*pi*f;

% air density, frequency etc
Constant = 1.224*(w^2)/(2*pi);

% number of boxes (should be 100 x 100 if panel is made using linspace)
SzX = size(A,2);
SzY = size(A,1);

% divide panel into small boxes with area LxLy/SzXSxY
Dx = linspace(-Lx/2,Lx/2,SzX);
Dy = linspace(-Ly/2,Ly/2,SzY);

% Panel area
dS = (Lx*Ly)/(SzX*SzY);

% wavenumber
k = w/343;

out = 0;

% perform rayleigh integral (discrete)
for i = 1:SzX
    
    for h = 1:SzY
        
        R = sqrt((Dx(i) - x)^2 + (Dy(h) - y)^2 + z^2);
        
        out = out + dS*Constant*A(h,i)*exp(-1i*k*R)/R;
        
    end
    
end