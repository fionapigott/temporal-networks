% Fiona Pigott
% February 24 2014
% MATLAB V2012b

% DOESN'T WORK. A work in progress...

% Calculate the clustering coefficient and the 
% weighted clustering coefficient to identify triangles.

% OUTPUT:
% ci: 
% ciweighted: Matrix with dim (# nodes) X (# time steps).
%   ciweighted(i,m) = clustering coeff at node i at time m
%   Formula from TareasV1.pdf #11
%   formula (LaTeX) c_i^w = \frac{1}{s_i(k_i-1)}\sum_{j,h}
%                           \frac{(w_{ij}+w_{ih})a_{ij}a_{ih}a_{jh}}{2}

% Sum of all of the interactions of node i per time step
% s(i,m) = sum of all interactions of node i (summed along cols i)
tic
s = squeeze(sum(data,1));

triangles = 0;
trianglesweighted = 0;
ci = zeros(numnodes,nummat);
ciweighted = ci;
for m = 1:nummat
    for ii = 1:numnodes % for all nodes
        if s(ii,m) ~=0 % if the node has interactions with other nodes
            % calculate the weight factor
            weightFactor = 1/(s(ii,m));
            tripletFactor = 1/(k(ii,m)-1);
            for jj = 1:ii-1 % unique combos of i & j, i =/= j
                    if unweighted(ii,jj,m)
                        for h = 1:jj-1 % unique combos of i,j &h, j=/=h
                            if unweighted(ii,h,m) && unweighted(jj,h,m)
                                trianglesweighted = trianglesweighted +...
                                    (data(ii,jj,m) + data(ii,h,m))/2;
                                triangles = triangles + 1;
                            end
                        end
                    end
            end
            if triangles ~= 0 
                ciweighted(ii,m) = trianglesweighted*...
                    weightFactor*tripletFactor;
                triangleweighted = 0;
                ci(ii,m) = triangles*tripletFactor;
                triangles = 0;
            end
        end
    end
end
 toc
 