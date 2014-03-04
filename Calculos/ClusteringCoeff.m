% Fiona Pigott
% February 24 2014
% MATLAB V2012b

% Calculate the clustering coefficient and the 
% weighted clustering coefficient to identify triangles.
% Formulea form "The architecture of complex weighted networks"
% by Barrat et al

% OUTPUT:
% ci: Matrix with dim (# nodes) X (# time steps).
%   ci(i,m) - clustering coeff at node i at time m
%   formula:  c_i^w = \frac{1}{k_i(k_i-1)}\sum_{j,h} a_{ij}a_{ih}a_{jh}
% ciweighted: Matrix with dim (# nodes) X (# time steps).
%   ciweighted(i,m) = clustering coeff at node i at time m
%   formula (LaTeX) c_i^w = \frac{1}{s_i(k_i-1)}\sum_{j,h}
%                           \frac{(w_{ij}+w_{ih})a_{ij}a_{ih}a_{jh}}{2}
% CFunction: un-weighted clustering coefficient C(k) (same properties
%            as CwFunction below).
% CwFunction: weighted clustering coefficient averaged over
%             all vertices with degree k, C(k). vector with length
%             max(k) where CwFunction(k) = avg Cw for all nodes with 
%             degree = k. Note that this function does NOT preserve time
%             information.
% ClusterTime: unweighted (same as ClusterwTime, but unweighted)
% ClusterwTime: weighted clustering coefficient averaged over
%         all vertices in a give time step m. C(t). vector with length
%         (# time steps) where ClusterwTime(m) = avg Cw for all nodes at 
%         time m.

% Sum of all of the interactions of node i per time step
% s(i,m) = sum of all interactions of node i (summed along cols i)
%tic

triangles = 0;
trianglesweighted = 0;
ci = zeros(numnodes,nummat);
ciweighted = ci;
for m = 1:nummat
    for ii = 1:numnodes % for all nodes
        if s(ii,m) ~=0 % if the node has interactions with other 
            % nodes calculate the weight factor
            weightFactor = 1/(s(ii,m));
            degreeFactor = 1/(k(ii,m));
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
                trianglesweighted = 0;
                ci(ii,m) = triangles*tripletFactor*degreeFactor;
                triangles = 0;
            end
        end
    end
end

civector = reshape(ci,numnodes*nummat,1);
kvector = reshape(k,numnodes*nummat,1);
CFunction = accumarray(kvector(kvector~=0),civector(kvector~=0))./...
    accumarray(kvector(kvector~=0),1);

ciweightedvector = reshape(ciweighted,numnodes*nummat,1);
CwFunction = accumarray(kvector(kvector~=0),ciweightedvector(kvector~=0))./...
    accumarray(kvector(kvector~=0),1);

ClusterTime = mean(ciweighted);
ClusterwTime = mean(ci);

%toc
clear triangles trianglesweighted weightFactor tripletFactor degreeFactor
clear s
 