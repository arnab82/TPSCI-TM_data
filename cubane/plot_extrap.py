import os
import sys
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
import matplotlib.lines as mlines
import scipy
from matplotlib.ticker import FuncFormatter
file = os.path.splitext(sys.argv[1])[0]

print(" Working with file: ", file)


np.set_printoptions(suppress=True, precision=6, linewidth=1500)

data_path = file+'.csv'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader)
    data = list(reader)

conversion = 219474.63 # cm-1
# conversion = 1 # hartrees
# conversion = 1000 # mH 
# conversion =27.2114 # eV
n_extrap_points = int((len(data[0]))//2)
print(" Number of extrapolation points: ", n_extrap_points)
energy_var = {}
energy_pt2 = {}
num_thresh = 6
 


for i in range(len(data) - 1): 
    energy_var[i] = np.array([float(a) * conversion for a in data[i+1][0:n_extrap_points]])
    energy_pt2[i] = np.array([float(a) * conversion for a in data[i+1][n_extrap_points:2*n_extrap_points]])


emin = 0
for i in energy_pt2:
    emin = min(emin, np.min(energy_pt2[i]))
for i in energy_pt2:
    energy_pt2[i] -= emin
    print(energy_pt2[i])
    energy_var[i] -= emin
    print(energy_var[i])

print(" Energy shift: %12.8f" %emin)

blue = '#3E6D9C'
orange = '#FD841F'
red_orange = '#E14D2A'
dark_blue = '#001253'

cb = ['#000000']
cb.extend([i for i in plt.rcParams['axes.prop_cycle'].by_key()['color']])
print(cb)
# cc = ['#FF0000', '#0000FF', '#FF00FF',  '#800000', '#008000', '#000080', '#808000', '#800080', '#008080', '#808080', '#C0C0C0', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#800000', '#008000', '#000080', '#808000', '#800080', '#008080', '#808080', '#C0C0C0','#FF1493', '#00FF7F', '#1E90FF', '#FF4500', '#9400D3', '#FF8C00', '#00CED1', '#FF69B4', '#00BFFF', '#FFD700', '#8A2BE2', '#32CD32', '#FF6347', '#00FFFF', '#FFA500', '#6A5ACD', '#FF00FF', '#7FFF00', '#9932CC', '#FF8C69', '#00FA9A', '#BA55D3', '#FF7F50', '#ADFF2F', '#8B008B', '#FF4500']
# '#00008B' ground state color
cc=['#00008B','#808000','#808000','#808000','#008080','#008080','#008080','#800080','#800080','#800080','#008000','#008000','#008000','#008000','#008000','#008000','#4B0082','#4B0082','#4B0082']
fig, ax = plt.subplots()

#this sets the number of decimal points on axis energies
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
extrap = []
xmax = -750
xmin = -800

ymax = 1500
ymin = -600



print(energy_pt2)
print(energy_var)
for key in energy_var:
    x = energy_pt2[key] - energy_var[key]
    z = energy_pt2[key]
    y = energy_var[key]

    m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    m2, b2, r_value, p_value, std_err = scipy.stats.linregress(x, z)
    
    extrap.append(b)

    plt.rcParams.update({'font.size': 10})

    ymin = min(ymin, b)
    ymax = max(ymax, m*xmin+b)
    print(x,y)
    print(x,z)
    # if key >=1:
    #     break
    for j in range(len(x)):
        ax.plot(x[j], y[j]+446.5539526023142, marker='o', linestyle='-', markersize=8, color=cc[key])
        ax.plot(x[j], z[j]+446.5539526023142, marker='x', linestyle=' ', markersize=8, color=cc[key])
    
    x2 = np.array([-1,0])*conversion
    line = m*x2+b+446.5539526023142
    print("y axis values are", line)
    print("b value is ", b)
    ax.plot(x2, line, alpha=1.0, color=cc[key], linestyle='-', linewidth=1.5  )
    line = m2*x2+b +446.5539526023142 
    print("y axis values are", line)                                  
    ax.plot(x2, line, alpha=0.5, color = cc[key], linestyle='--', linewidth=1.5)
    print("Extrapolated Result: %14.8f"% ((b+emin)/conversion))
    print("R^2                : %14.8f"% r_value)
    print("Var root",key,y)
    print("PT  root",key,z)
    #print("DIFFF",x)
    #print("color",cb[key])
var_marker = mlines.Line2D([], [], color='grey', marker='o', linestyle='None',
                          markersize=10, label='Variational')
pt2_marker = mlines.Line2D([], [], color='grey', marker='x', linestyle='None',
                          markersize=10, label='PT2')

print("x: ", (xmin, xmax))
print("y: ", (ymin, ymax))
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)
# Set y-axis ticks and labels
ax.set_xlabel('$\Delta$E$_{PT2}$ (eV) ')
ax.set_ylabel('Shifted Energy (eV) ')
# ax.set_xlabel('$\Delta$E$_{PT2}$ (mH) ')
# ax.set_ylabel('Shifted Energy (mH) ')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

formatter = FuncFormatter(lambda x, _: f'{x:.2f}')
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)

#ax.set_yticklabels([])
# ax.yaxis.set_label_position("right")
# ax.yaxis.tick_right()
#ax.legend()
# ax.set_title('Absolute Energy, $E(S)$ (mH) ')
#ax.set_title("Tetracence Tetramere \n (40o, 40e) \n ε = {0.0005, 0.0007, 0.001, 0.005, 0.01} \n 31 roots")
#do some legend stuff or comment out
# black_patch = mpatches.Patch(color=dark_blue, label='Ground')
# blue_patch = mpatches.Patch(color='#008000', label='Triplets')
# green_patch = mpatches.Patch(color='#4B0082', label='Singlets')
# red_patch = mpatches.Patch(color='#800080', label='Biexcitons')
# ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='center right',bbox_to_anchor=(1, 0.6))
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='upper center')
# ax.legend()
# handles, labels = ax.get_legend_handles_labels()
# filtered_handles = [handle for handle in handles if handle.get_linestyle() == '-']
# filtered_labels = [label for label, handle in zip(labels, handles) if handle.get_linestyle() == '-']
# ax.legend(filtered_handles, filtered_labels, loc='upper left',bbox_to_anchor=(1, 1))
# var_marker = mlines.Line2D([], [], color='grey', marker='o', linestyle='None',
#                           markersize=10, label='Variational')
# pt2_marker = mlines.Line2D([], [], color='grey', marker='x', linestyle='None',
#                           markersize=10, label='PT2')
# fig = plt.gcf()
fig.set_size_inches(4,7)
fig.savefig(file+'_extrap2.pdf', dpi=300, bbox_inches='tight')

plt.show()

print(extrap)
