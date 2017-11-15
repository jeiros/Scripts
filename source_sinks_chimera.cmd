colordef transgreen 0 1 0 .3
colordef transred 1 0 0 .3
colordef transblue 0 0 1 .3
preset apply pub 1
sel :65,66,67,69,71,73,75,76,280,420
~disp
disp sel
cofr :420
color transred #0-4
color transblue #5-9

~ribbon
ribbon :1-161
ribbon :249-320
cofr :420



# for morphng
~disp
sel :65,66,67,69,71,73,75,76,280,420,271,272
disp sel
color transblue :1-161; color transred :249-320
