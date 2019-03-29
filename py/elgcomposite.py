import jplus
import numpy as np
import matplotlib.pyplot as plt
jplus.plotting.setup_text_plots(fontsize=14,usetex=True)
spec_elg = jplus.datasets.fetch_eboss_elg_composite()


z = 0.76
#spec_elg['w'] *= (1 + z)

ewarr = np.linspace(10,100,num=5)
color = plt.cm.Paired(np.linspace(0.1,0.9,5))
k = 0
for i in ewarr:
  specline = jplus.datasets.add_line_to_spec(spec_elg,EW=i,obsframe=True, z=z)
  plt.plot(specline['w'],specline['flux'][:,0],color=color[k])
  print 'EW = ',i
  k += 1

plt.xlim([3000,9000])    
plt.show()  
