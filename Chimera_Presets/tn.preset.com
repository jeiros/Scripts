~disp
color blue :1-161
color green :162-248
color red :249-419
disp :CAL
repr sphere :CAL
color silver :CAL
repr sphere ~protein
disp ~protein
