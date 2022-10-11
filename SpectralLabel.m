%% labels = SpectralLabel(X,gamma)
% Creates bit labels from Fiedler vector loosely based on spectral
% clustering
function  labels = SpectralLabel(X,gamma)
% Assumes Gaussian, i.e., squared euclidean distance
EucD = pdist2(X,X,'squaredeuclidean');

% affinity matrix
A = exp(-gamma*EucD) - eye(size(X,1));
% diagionals
D = sum(A,1);

%lagrangian
L = diag(D) - A; 
%normalisation
L = 1./sqrt(D).' .* L .* 1./sqrt(D); 

% get eigenvectors
[x, ~] = eig(L);
% undo normalisation
vectors = 1./sqrt(D).' .* x;
% get order of the Fiedler vector
[~,idx] = sort(vectors(:,2));

% % deduct 1 to got from 1-indexed to binary label
% labels = idx-1;
labels(idx,1) = 0:size(X,1)-1;
end