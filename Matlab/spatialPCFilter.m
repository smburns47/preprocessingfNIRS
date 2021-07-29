function filtered = spatialPCFilter(d, coords)
%adapted from code by Zhang, Noah, & Hirsch, 2016
%needs wide coverage area to work appropriately, so don't use if less than 
%~20 channels
%TO DO: simulations to test exactly how much coverage needed

data = d;
coords = [coords; coords];
[U,S,V] = svd(data);
[th,phi,~] = cart2sph(coords(:,1),coords(:,2),coords(:,3));
thphi = nan(size(coords,1),2);
thphi(:,1) = th;
thphi(:,2) = phi;
distmatrix = squareform(pdist(thphi,@arclen));
kernel = exp(- distmatrix.^2 / 50);
for y=1:size(kernel,2)
    kernel(:,y) = kernel(:,y) / sum(kernel(:,y));
end
Vsmooth = nan(size(V));
for pc = 1:size(coords,1)
    for ch = 1:size(coords,1)
        Vsmooth(ch,pc) = V(:,pc)' * kernel(:,ch);
    end
end
Hglobal = U*S*Vsmooth';
filtered = data - Hglobal;
end

function l=arclen(XI,XJ)
s=referenceSphere('Unit Sphere');
l=distance('gc',XI*180/pi,XJ*180/pi,s);
end