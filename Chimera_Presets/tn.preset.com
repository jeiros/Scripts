~disp
colordef transwhite 1 1 1 .3
colordef transred 1 0 0 .1
colordef transgreen 0 1 0 .1
colordef transblue 0 0 1 .1
color blue :1-161
color green :162-248
color red :249-419
disp :CAL
repr sphere :CAL
color silver :CAL
repr sphere ~protein
disp ~protein
~sel
