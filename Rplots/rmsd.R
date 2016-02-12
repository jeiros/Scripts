

r1 <- read.table("rmsd_WT_1.000-750ns.dat_time", header = T)
r2 <- read.table("rmsd_WT_2.000-750ns.dat_time", header = T)
r3 <- read.table("rmsd_WT_3.000-750ns.dat_time", header = T)
r4 <- read.table("rmsd_WT_4.000-750ns.dat_time", header = T)
r4 <- read.table("rmsd_WT_5.000-750ns.dat_time", header = T)




ggplot() + theme_bw(15) +
  geom_line(data=r1, aes(x=X0, y=rmsd, color = 'Run1'), size = 1) +
  geom_line(data=r2, aes(x=X0, y=rmsd, color = 'Run2'), size = 1) +
  geom_line(data=r3, aes(x=X0, y=rmsd, color = 'Run3'), size = 1) +
  geom_line(data=r4, aes(x=X0, y=rmsd, color = 'Run4'), size = 1) +
  labs(y="RMSD (Å)", x="Time (ns)") +
  scale_x_continuous(breaks = c(0,150,300,450,600,750)) +
  scale_color_discrete("WT run") +
  theme(legend.justification=c(1,0), legend.position=c(.95, .05)) +
  theme(legend.background = element_rect(fill='white', colour='black'))