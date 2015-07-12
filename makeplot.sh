#!/bin/bash

for i in `ls TnCI_*.dat`; do
	matlab -nodisplay <<- EOF
		A = dlmread('$i','\t', [1 0 -1 2]); 
		jmin = min(A(:,1));
		kmin = min(A(:,2));
		length = size(A,1);

		for i = 1:length
		    j = A(i,1);
		    k = A(i,2);
		    M((j-jmin+1),(k-kmin+1)) = A(i,3);
		end
		imagesc(M)
		xlabel('TnI')
		ylabel('TnC')
		colorbar
		set(gca,'YDir','normal')
		print('$i.png','-dpng')
		exit
	EOF
done