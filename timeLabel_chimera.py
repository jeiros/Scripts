# Issue this command inside Chimera
# 2dlabels create timer text "0 ns" color black ypos .9 xpos .9 size 36

# Issue this per-frame script
from chimera import runCommand
runCommand("2dlabels change timer text '%.1f ns'" % (mdInfo['frame']*0.02))
