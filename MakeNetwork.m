// Load the mapping file (msms.mapping_)  and the Trans matrix and the labels file



SortMatrix = zeros(128);
veccount = zeros(128,1);
// each line is a cluster that contains all the frames in that cluster
// and the first line contains the number of frames in that cluster
for i = 1:NoFrames
   SortMatrix(Labels_file(1,i)+1, veccount(Labels_file(1,i)+1, 1) +1) = i;
   veccount(Labels_file(1,i)+1, 1) = veccount(Labels_file(1,i)+1, 1) +1;
end

count = 1;
for i = 1:Dim(Transmat)
for j = 1:Dim(Transmat)
if (Transmat_file(i,j) ~= 0 && i~=j)
network(count,1) = SortMatrix(Mapping_file(i,1) + 1,1);
network(count,2) = Transmat_file(i,j);
network(count,3) = SortMatrix(Mapping_file(j,1) + 1,1);
network(count,4) = veccount(Mapping_file(i,1) + 1,1);
network(count,5) = veccount(Mapping_file(j,1) + 1,1);
// network file has 4 lines 
// 1st line is the 1st fr
count = count +1;
end
end
end

dlmwrite('OutputName', network, 'delimiter', '\t', 'precision', 8)
