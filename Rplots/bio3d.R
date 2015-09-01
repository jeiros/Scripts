library(bio3d)
library(ggplot2)

setwd("/home/je714/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs/CTnI_hmr/run1/S1P")



pdb <- read.pdb("output.pdb")
traj_files <- list.files(".", "05_Prod")


trajs <- read.ncdf(traj_files, stride = 1)

ca.inds <- atom.select(pdb, elety="CA") 
# ca.inds is a list containing atom 
# and xyz numeric indices that we can 
# now use to superpose all frames of the trajectory 
# on the selected indices (in this
# case corresponding to all alpha Carbon atoms)


xyz <- fit.xyz(fixed = pdb$xyz, mobile = trajs, 
			   fixed.inds = ca.inds$xyz, 
			   mobile.inds = ca.inds$xyz)
# performs the actual superposition and stores the new coordinates
# in the matrix object xyz. Note that the dimensions (i.e. number 
# of rows and columns, which correspond to frames and coordinates
# respectively) of xyz match those of the input trajectory

rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frames <- dim(trajs)[1]
time <- 50 * length(traj_files) # Each file has 50 ns of simulation
rd_dataframe <- data.frame(rd, seq(1/(frames/time), time, 1/(frames/time)))
colnames(rd_dataframe) <- c("RMSD", "Time")
p <- ggplot() + geom_line(data = rd_dataframe, aes(x = Time, y = RMSD), colour = "black") + 
	ylab("RMSD (A)") + xlab("Time (ns)")
print(p)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
time1 <- seq(0.0125, 50, 0.0125)
time2 <- seq(50.02, 200, 0.02)
rd_dataframe1 <- data.frame(rd[1:4000], time1)
rd_dataframe2 <- data.frame(rd[4001:11500], time2)
colnames(rd_dataframe1) <- colnames(rd_dataframe2) <- c("RMSD", "Time")
p <- ggplot() + geom_line(data = rd_dataframe1, aes (x = Time, y = RMSD), 
	colour = "black") + geom_line(data = rd_dataframe2, aes(x = Time, 
		y = RMSD), colour = "black") + ylab("RMSD (Å)") + xlab("Time (ns)")
print(p)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#ggsave("ggsave2.png", width = 8.5, unit = "cm", dpi = 300)

#RMSD histogram
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)

#RMSF plot
rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")

#PCA
pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)))


# Below we perform a quick clustering in PC-space to further highlight these distinct conformers.
hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)

# examine the contribution of each residue to the first two principal components.

plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")


# 
p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="pc2.pdb")
