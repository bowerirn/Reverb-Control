%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code calculates modal resonant frequencies up to Mx by My for a 
%   plate of size Lx by Ly.The material properties are included in this 
%   code, as well as the plate thickness h.
%
%   ALUMINUM         G. GLASS          STEEL            ACRYLIC
%   E = 68.9e9;      E = 71.5e9;       E = 200e9;       E = 3.2e9
%   v = 0.334;       v = 0.21          v = 0.291        v = 0.35
%   rho = 2700;      rho = 2450        rho = 7850       rho = 1180 
%
%   DAA edited 6/3/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = single_clamped_mode_freqs(Lx,Ly,h,E,v,rho,Mx,My)



A = E*h^2;
B = 12*rho*(1-v^2);
D = (pi^2)*sqrt(A/B);

Lx2 = Lx^2;
Ly2 = Ly^2;

% change to 1:Mx 1:My to make array of res freqs
for m = Mx
    for n = My
        
        dm = ((((n*Lx)/(m*Ly))^2) + 2)^(-1);
        dn = ((((m*Ly)/(n*Lx))^2) + 2)^(-1);
        
        C = ((m + dm)^2)/Lx2 + ((n + dn)^2)/Ly2;
        % Change to M(m,n)
        M = D*C;
        
    end
end