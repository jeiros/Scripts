~disp
colordef transwhite 1 1 1 .4
colordef transblack 1 1 1 .4
colordef transred 1 0 0 .3
colordef transgreen 0 1 0 .3
colordef transblue 0 0 1 .3
color blue :1-161
color green :162-248
color red :249-419
disp :CAL
repr sphere :CAL
color silver :CAL
repr sphere ~protein
disp ~protein
~sel
