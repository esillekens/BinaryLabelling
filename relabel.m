%% function X = relabel(X,map)
% applies the natural bit mapping to constellation X 
% if map is the bitlabel
function X = relabel(X,label)
    X(label+1,:) = X;
end