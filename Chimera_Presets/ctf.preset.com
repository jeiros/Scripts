~display
background solid black
colordef transwhite 1 1 1 .3
# F-actin
ribbon :1-6016
color transwhite :1-6016

# Tropomyosin coil n1
ribbon :6033-6712
color purple :6033-6712

# Tropomyosin coil n2
ribbon :7336-8015
color purple :7336-8015


# Troponin complex 1
# TnT
ribbon :6713-6961
color green :6713-6961
# TnI
ribbon :6962-7171
color red :6962-7171
# TnC
ribbon :7172-7332
color blue :7172-7332

# Troponin complex 2
# TnT
ribbon :8016-8264
color green :8016-8264
# TnI
ribbon :8265-8474
color red :8265-8474
# TnC
ribbon :8475-8635
color blue :8475-8635


# CALCIUMS
disp :CAL
repr sphere :CAL
color silver :CAL

# ADPs
disp :ADP
repr sphere :ADP
