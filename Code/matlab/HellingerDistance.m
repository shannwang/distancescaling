function H=HellingerDistance(p,q,de)
% calculate the Hellinger distance
% p and q are two distributions
% de is the width of bins
H=sqrt(1-sqrt(p)*sqrt(q')*de);
end