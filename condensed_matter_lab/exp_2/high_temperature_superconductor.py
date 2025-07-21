import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data = np.loadtxt("ybco_down (7).all")

temperature = data[:,0]
magnetic_field = data[:,1]
area = data[:,2]
rho = data[:,3]
current = data[:,4]
voltage = data[:,5]
resistance = data[:,6]
u_offset = data[:,7]
time = data[:,8]

i = 0
positive_resistance = []
positive_temperature = []
negative_resistance = []
negative_temperature = []
avg_negative_temperature = []
avg_negative_resistance = []
avg_positive_temperature = []
avg_positive_resistance = []

print(temperature[0], temperature[-1])


# raw curve
fig, ax0 = plt.subplots()
ax0.plot(temperature, [val*1000 for val in resistance], linewidth = 1, c="black")
# Set the x-axis and y-axis labels for ax1
ax0.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax0.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax0.set_title("Heat Down", size=16)
# margin: left,bottom,right,top
ax0.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax0.tick_params(axis='x', labelsize=20)
ax0.tick_params(axis='y', labelsize=20)
ax0.set_ylim(-.2,2.8)
ax0.grid(True, linestyle='--', alpha=0.7)
#ax1.legend()
#plt.show()
ax0.figure.savefig("raw_resistance_down.svg", format = "svg")





while i < len(current)-1:
    if current[i]<0:
        j = i
        count1=0
        sum_neg_resistance = 0
        sum_neg_temperature = 0
        while current[j]<0 and j < len(current)-1:
            count1 += 1
            sum_neg_resistance += resistance[j]
            sum_neg_temperature += temperature[j]
            j += 1
        i = j
        #print((i-(int(count1/2)-1))-(i+(int(count1/2))-1))
        # average between positive and negative points
        avg_negative_temperature.append(np.mean(temperature[i-(int(count1/2)-1):i+(int(count1/2))-1]))
        avg_negative_resistance.append(np.mean(resistance[i-(int(count1/2)-1):i+(int(count1/2))-1]))
        # average of negatives points
        negative_resistance.append(sum_neg_resistance/count1)
        negative_temperature.append(sum_neg_temperature/count1)

    elif current[i] >= 0:
        k = i
        count2 = 0
        sum_pos_resistance = 0
        sum_pos_temperature = 0
        while current[k]>=0  and k < len(current)-1:
            count2 += 1
            sum_pos_resistance += resistance[k]
            sum_pos_temperature += temperature[k]
            k += 1
        i = k
        #print(i-(int(count2/2)-1),(i+(int(count2/2))))
        # average between positive and negative points
        avg_positive_temperature.append(np.mean(temperature[i-(int(count2/2)-1):i+(int(count2/2))]))
        avg_positive_resistance.append(np.mean(resistance[i-(int(count2/2)-1):i+(int(count2/2))]))
        # average of positive points
        positive_resistance.append(sum_pos_resistance / count2)
        positive_temperature.append(sum_pos_temperature/count2)

"""
plt.figure(figsize=(8, 6))
#plt.scatter(positive_temperature,positive_resistance, marker = "+", c="red")
#plt.scatter(negative_temperature,negative_resistance, marker = "+")
#plt.scatter(avg_negative_temperature,avg_negative_resistance, marker = "+", c="black")
#plt.scatter(avg_positive_temperature,avg_positive_resistance, marker = "+", c="b")
#plt.show()

# whole curve
fig, ax1 = plt.subplots()
ax1.scatter(avg_negative_temperature,[val*1000 for val in avg_negative_resistance], marker = ".", c="black")
ax1.scatter(avg_positive_temperature,[val*1000 for val in avg_positive_resistance], marker = ".", c="black")
# Set the x-axis and y-axis labels for ax1
ax1.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax1.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax1.set_title("Cool Down", size=16)
# margin: left,bottom,right,top
ax1.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
#ax1.legend()
#plt.show()
ax1.figure.savefig("resistance_down.svg", format = "svg")

# zoomed curve
fig, ax2 = plt.subplots()
ax2.scatter(avg_negative_temperature,[val*1000 for val in avg_negative_resistance], marker = ".", c="black")
ax2.scatter(avg_positive_temperature,[val*1000 for val in avg_positive_resistance], marker = ".", c="black")
# Set the x-axis and y-axis labels for ax1
ax2.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax2.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax2.set_title("Cool Down", size=16)
# margin: left,bottom,right,top
ax2.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(72,150)
ax2.set_ylim(-.2,1.5)
#ax1.legend()
#plt.show()
ax2.figure.savefig("zoomed_resistance_down.svg", format = "svg")
"""

"""

# Create the inset plot of resistivity
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot the whole curve in the main plot
ax1.scatter(avg_negative_temperature, [val * 1000 for val in avg_negative_resistance], marker=".", c="black")
ax1.scatter(avg_positive_temperature, [val * 1000 for val in avg_positive_resistance], marker=".", c="black")
ax1.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax1.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax1.set_title("Cool Down", size=16)
# margin: left,bottom,right,top
ax1.set_position([0.12, .13, 0.82, 0.82])
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)

# Create the inset plot with more space at the bottom
inset_bbox = [0.12, 0.2, 0.82, 0.82]  # Adjust the values as needed
axins = inset_axes(ax1, width="40%", height="40%", loc='lower right', bbox_to_anchor=inset_bbox, bbox_transform=fig.transFigure)

# Plot the zoomed curve in the inset plot
axins.scatter(avg_negative_temperature, [val * 1000 for val in avg_negative_resistance], marker=".", c="black")
axins.scatter(avg_positive_temperature, [val * 1000 for val in avg_positive_resistance], marker=".", c="black")
axins.set_xlim(72, 150)
axins.set_ylim(-0.2, 1.5)
axins.set_position([0.12, .2, 0.82, 0.82])

# Set ticks for the inset plot
axins.set_xticks([80, 100, 120, 140])
axins.set_yticks([0, 0.5, 1, 1.5])

# Add a rectangle indicating the zoomed area in the main plot
ax1.indicate_inset_zoom(axins)

# Save the figure
plt.savefig("resistance_down_inset.svg", format="svg")






# resistivity curve
fig, ax3 = plt.subplots()
ax3.scatter(avg_negative_temperature,[val*1000000*2*np.pi*.0026 for val in avg_negative_resistance], marker = ".", c="black")
ax3.scatter(avg_positive_temperature,[val*1000000*2*np.pi*.0026 for val in avg_positive_resistance], marker = ".", c="black")
# Set the x-axis and y-axis labels for ax1
ax3.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax3.set_ylabel('Resistivity (µΩ.m)', size=20, labelpad=10)
ax3.set_title("Cool Down", size=16)
# margin: left,bottom,right,top
ax3.set_position([0.15, .17, 0.84, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(72,150)
ax3.set_ylim(-3,25)
#ax1.legend()
#plt.show()
ax3.figure.savefig("resistivity_down.svg", format = "svg")
"""




# midpoint
def midpoint(x1, y1, x2, y2):
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    return mid_x, mid_y

# mid point curve
fig, ax4 = plt.subplots()
merged_temperature = np.vstack((avg_negative_temperature, avg_positive_temperature)).ravel(order='F')
merged_resistance = np.vstack((avg_negative_resistance, avg_positive_resistance)).ravel(order='F')
point1 = (merged_temperature[-6],merged_resistance[-6]*1000)
point2 = (merged_temperature[-9],merged_resistance[-9]*1000)
midpoint_result = midpoint(point1[0], point1[1], point2[0], point2[1])
T_c = "{:.2f}".format(midpoint_result[0])
print("start:",merged_temperature[-6], "end:", merged_temperature[-9], "T_c:",T_c)
# main data
ax4.plot(merged_temperature, [val*1000 for val in merged_resistance], marker =".", linestyle="-", c="black")
# temperature
ax4.axvline(merged_temperature[-6], color='red', linestyle='--', alpha = .5, label=f"T = {point1[0]:.2f}K ,{point2[0]:.2f}K")
ax4.axvline(merged_temperature[-9], color='red', linestyle='--', alpha = .5)
# resistance
ax4.hlines(merged_resistance[-6]*1000, xmin=70, xmax=125, color='b', linestyle='--', alpha = .5, label=f"T = {point1[1]:.2f}mΩ ,{point2[1]:.2f}mΩ")
ax4.hlines(merged_resistance[-9]*1000, xmin=70, xmax=125, color='b', linestyle='--', alpha = .5)
# mid point
ax4.axvline(midpoint_result[0], color='g', linestyle='--', alpha = .5)
ax4.hlines(midpoint_result[1], xmin=70, xmax=125, color='g', linestyle='--', alpha = .5, label=f"mid point ($T_c$ = {T_c}K)")
# Set the x-axis and y-axis labels for ax1
ax4.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax4.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
# margin: left,bottom,right,top
ax4.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
ax4.set_title("Cool Down", size=16)
ax4.set_xlim(65,150)
ax4.set_ylim(-.2,1.5)
ax4.legend()
ax4.grid(True, linestyle='--', alpha=0.7)
#plt.show()
ax4.figure.savefig("mid_point_down.svg", format = "svg")






# Create the inset plot of resistance
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot the whole curve in the main plot
ax1.plot(merged_temperature, [val * 1000 for val in merged_resistance], marker=".",linestyle="-", c="black")
ax1.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax1.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax1.set_title("Cool Down", size=16)
# margin: left,bottom,right,top
ax1.set_position([0.12, .13, 0.82, 0.82])
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.grid(True, linestyle='--', alpha=0.7)

# Create the inset plot with more space at the bottom
inset_bbox = [0.12, 0.2, 0.82, 0.82]  # Adjust the values as needed
axins = inset_axes(ax1, width="40%", height="40%", loc='lower right', bbox_to_anchor=inset_bbox, bbox_transform=fig.transFigure)

# Plot the zoomed curve in the inset plot
axins.plot(merged_temperature, [val * 1000 for val in merged_resistance], marker=".",linestyle="-", c="black")
axins.set_xlim(72, 150)
axins.set_ylim(-0.2, 1.5)
axins.set_position([0.12, .2, 0.82, 0.82])

# Set ticks for the inset plot
axins.set_xticks([80, 100, 120, 140])
axins.set_yticks([0, 0.5, 1, 1.5])
axins.grid(True, linestyle='--', alpha=0.7)

# Add a rectangle indicating the zoomed area in the main plot
ax1.indicate_inset_zoom(axins)

# Save the figure
plt.savefig("resistance_down_inset.svg", format="svg")






# resistivity curve
fig, ax3 = plt.subplots()
ax3.plot(merged_temperature,[val*1000000*2*np.pi*.0026 for val in merged_resistance], marker=".",linestyle="-", c="black")
# Set the x-axis and y-axis labels for ax1
ax3.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax3.set_ylabel('Resistivity (µΩ.m)', size=20, labelpad=10)
ax3.set_title("Cool Down", size=16)
# margin: left,bottom,right,top
ax3.set_position([0.15, .17, 0.84, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(72,150)
ax3.set_ylim(-3,25)
ax3.grid(True, linestyle='--', alpha=0.7)
#ax1.legend()
ax3.figure.savefig("resistivity_down.svg", format = "svg")







# ------------------------ xrd ---------------------

bfxrd_data = np.loadtxt("bs-1.txt")
afxrd_data = np.loadtxt("as_1.txt")

#print(bfxrd_data[:,0])
# xrd curve
fig, ax5 = plt.subplots()

#print(len(afxrd_data[:,0]),len(afxrd_data[:,1]))
ax5.plot(bfxrd_data[:,0],bfxrd_data[:,1], linewidth = .5, c = "black")
ax5.plot(afxrd_data[:,0],[value+800 for value in afxrd_data[:,1]], linewidth = .5, c= "r")
# before centering
ax5.text(19.7,400 , "+(110)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(26.5,800 , "+(111)", ha='center', va='center', fontsize=9, rotation = 0, c = "b")
ax5.text(25.5,500 , "+(021)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(27.7,380 , "+(002)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(29.6,570 , "#(111)", ha='center', va='center', fontsize=9, rotation = 90, c = "r")
ax5.text(33.6,500 , "+(112)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(34.9,700 , "+(130)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(36.3,480 , "*(-111)", ha='center', va='center', fontsize=9, rotation = 90, c = "g")
ax5.text(39,500 , "*   (-111)", ha='center', va='center', fontsize=9, rotation = 90, c = "g")
ax5.text(40,380 , "+(220)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(42.2,700 , "+(221)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(43.3,390 , "+(041)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(45,500 , "+(132)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(47,400 , "+(113)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(49,500 , "#(220)", ha='center', va='center', fontsize=9, rotation = 90, c = "r")
ax5.text(55.8,400 , "#(311)", ha='center', va='center', fontsize=9, rotation = 90, c = "r")
ax5.text(57.8,400 , "+(151)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(66.2,400 , "*(-311)", ha='center', va='center', fontsize=9, rotation = 90, c = "g")
ax5.text(68.3,400 , "+(332)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
ax5.text(71.2,400 , "+(400)", ha='center', va='center', fontsize=9, rotation = 90, c = "b")
# after centering
ax5.text(23,1350 , "x (010)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(27.8,1300 , "x (102)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(32,2200 , "x (103)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(35.3,2950 , "x (110)", ha='center', va='center', fontsize=9, rotation = 0, c = "c")
ax5.text(36.4,1250 , "x (112)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(38.5,1400 , "x (005)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(40.4,1500 , "x (113)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(46.7,1550 , "x (020)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(48,1350 , "x (200)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(51.6,1250 , "x (115)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(58.3,1750 , "x (123)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(59.2,1350 , "x (213)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(68.8,1350 , "x (220)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")
ax5.text(77.8,1300 , "x (130)", ha='center', va='center', fontsize=9, rotation = 90, c = "c")



ax5.add_patch(Rectangle((62, 2100), 16, 800,
             edgecolor = 'black',
             facecolor = 'none',
             fill=False,
             lw=.5))

# keys
ax5.text(65,2750, "+ BaCO$_3$", fontsize=12, c="b")
ax5.text(65,2550, "# Y$_2$O$_3$", fontsize=12, c="r")
ax5.text(65,2350, "* CuO", fontsize=12, c="g")
ax5.text(65,2150, "x YBCO", fontsize=12, c="c")



# Set the x-axis and y-axis labels for ax1
ax5.set_xlabel('2$\\theta$ (degree)', size=16, labelpad=10)
ax5.set_ylabel('Intensity (arb. units)', size=16, labelpad=10)
# margin: left,bottom,right,top
ax5.set_position([0.16, .15, 0.815, 0.83])
# Set the font size for x-ticks and y-ticks for ax1
ax5.tick_params(axis='x', labelsize=16)
ax5.tick_params(axis='y', labelsize=16)
ax5.set_xlim(10,80)
#ax5.set_ylim(-.2,1.5)
ax5.xaxis.set_ticks(np.arange(10, 90, 10))
#ax5.legend()
#ax5.legend(handles=[line1, line2, line3, line4])
#plt.show()
ax5.figure.savefig("xrd.svg", format = "svg")



