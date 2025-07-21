import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

data = np.loadtxt("ybco_up (2).all")

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




# raw curve
fig, ax0 = plt.subplots()
ax0.plot(temperature, [val*1000 for val in resistance], linewidth = 1, c="black")
# Set the x-axis and y-axis labels for ax1
ax0.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax0.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax0.set_title("Heat Up", size=16)
# margin: left,bottom,right,top
ax0.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax0.tick_params(axis='x', labelsize=20)
ax0.tick_params(axis='y', labelsize=20)
ax0.set_ylim(-.2,2.8)
ax0.grid(True, linestyle='--', alpha=0.7)
#ax1.legend()
#plt.show()
ax0.figure.savefig("raw_resistance_up.svg", format = "svg")



print(temperature[0], temperature[-1])


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
        # average between positive and negative points
        avg_positive_temperature.append(np.mean(temperature[i-(int(count2/2)-1):i+(int(count2/2))]))
        avg_positive_resistance.append(np.mean(resistance[i-(int(count2/2)-1):i+(int(count2/2))]))
        # average of positive points
        positive_resistance.append(sum_pos_resistance / count2)
        positive_temperature.append(sum_pos_temperature/count2)

plt.figure(figsize=(8, 6))
#plt.scatter(positive_temperature,positive_resistance, marker = "+", c="red")
#plt.scatter(negative_temperature,negative_resistance, marker = "+")
#plt.scatter(avg_negative_temperature,avg_negative_resistance, marker = "+", c="black")
#plt.scatter(avg_positive_temperature,avg_positive_resistance, marker = "+", c="b")
#plt.show()

"""
# whole curve
fig, ax1 = plt.subplots()
ax1.scatter(avg_negative_temperature,[val*1000 for val in avg_negative_resistance], marker = ".", c="black")
ax1.scatter(avg_positive_temperature,[val*1000 for val in avg_positive_resistance], marker = ".", c="black")
# Set the x-axis and y-axis labels for ax1
ax1.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax1.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax1.set_title("Heat Up", size=16)
# margin: left,bottom,right,top
ax1.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
#ax1.legend()
#plt.show()
ax1.figure.savefig("resistance_up.svg", format = "svg")


# zoomed curve
fig, ax2 = plt.subplots()
ax2.scatter(avg_negative_temperature,[val*1000 for val in avg_negative_resistance], marker = ".", c="black")
ax2.scatter(avg_positive_temperature,[val*1000 for val in avg_positive_resistance], marker = ".", c="black")
# Set the x-axis and y-axis labels for ax1
ax2.set_title("Heat Up", size=16)
ax2.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax2.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
# margin: left,bottom,right,top
ax2.set_position([0.16, .17, 0.82, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(72,150)
ax2.set_ylim(-.2,1.5)
ax2.legend()
#plt.show()
ax2.figure.savefig("zoomed_resistance_up.svg", format = "svg")
"""



"""
# Create the inset plot of resistivity
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot the whole curve in the main plot
ax1.scatter(avg_negative_temperature, [val * 1000 for val in avg_negative_resistance], marker=".", c="black")
ax1.scatter(avg_positive_temperature, [val * 1000 for val in avg_positive_resistance], marker=".", c="black")
ax1.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax1.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax1.set_title("Heat Up", size=16)
# margin: left,bottom,right,top
ax1.set_position([0.12, .13, 0.82, 0.82])
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.grid(True, linestyle='--', alpha=0.7)

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
axins.grid(True, linestyle='--', alpha=0.7)

# Add a rectangle indicating the zoomed area in the main plot
ax1.indicate_inset_zoom(axins)

# Save the figure
plt.savefig("resistance_up_inset.svg", format="svg")
"""







"""
# resistivity curve
fig, ax3 = plt.subplots()
#ax3.scatter(avg_negative_temperature,[val*1000000*2*np.pi*.0026 for val in avg_negative_resistance], marker = ".", c="black")
#ax3.scatter(avg_positive_temperature,[val*1000000*2*np.pi*.0026 for val in avg_positive_resistance], marker = ".", c="black")

# Use plot instead of scatter for line plot
ax3.plot(avg_negative_temperature, [val * 1000000 * 2 * np.pi * 0.0026 for val in avg_negative_resistance], marker=".", c="black", label='Negative')
ax3.plot(avg_positive_temperature, [val * 1000000 * 2 * np.pi * 0.0026 for val in avg_positive_resistance], marker=".", c="black", label='Positive')

# Set the x-axis and y-axis labels for ax1
ax3.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax3.set_ylabel('Resistivity (µΩ.m)', size=20, labelpad=10)
ax3.set_title("Heat Up", size=16)
# margin: left,bottom,right,top
ax3.set_position([0.15, .17, 0.84, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(72,150)
ax3.set_ylim(-3,25)
#ax1.legend()
#plt.show()
ax3.figure.savefig("resistivity_up_with_line.svg", format = "svg")
"""




# midpoint
def midpoint(x1, y1, x2, y2):
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    return mid_x, mid_y


# mid point curve
fig, ax4 = plt.subplots()

merged_temperature = np.vstack((avg_positive_temperature[:-1], avg_negative_temperature)).ravel(order='F')
merged_resistance = np.vstack((avg_positive_resistance[:-1], avg_negative_resistance)).ravel(order='F')
point1 = (merged_temperature[4],merged_resistance[4]*1000)
point2 = (merged_temperature[7],merged_resistance[7]*1000)
midpoint_result = midpoint(point1[0], point1[1], point2[0], point2[1])
T_c = "{:.2f}".format(midpoint_result[0])
print("start:",merged_temperature[4], "end:", merged_temperature[7], "T_c:",T_c)

# main data
ax4.plot(merged_temperature, [val*1000 for val in merged_resistance], marker=".", linestyle='-', c="black")
# temperature
ax4.axvline(merged_temperature[4], color='red', linestyle='--', alpha = .5, label=f"T = {point1[0]:.2f}K ,{point2[0]:.2f}K")
ax4.axvline(merged_temperature[7], color='red', linestyle='--', alpha = .5)
# resistance
ax4.hlines(merged_resistance[4]*1000, xmin=70, xmax=125, color='b', linestyle='--', alpha = .5, label=f"T = {point1[1]:.2f}mΩ ,{point2[1]:.2f}mΩ")
ax4.hlines(merged_resistance[7]*1000, xmin=70, xmax=125, color='b', linestyle='--', alpha = .5)
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
ax4.set_title("Heat Up", size=16)
ax4.set_xlim(65,150)
ax4.set_ylim(-.2,1.5)
ax4.grid(True, linestyle='--', alpha=0.7)
ax4.legend(loc='lower right')
#plt.show()
ax4.figure.savefig("mid_point_up.svg", format = "svg")



# resistivity curve
fig, ax3 = plt.subplots()

ax3.plot(merged_temperature, [val*1000000 * 2 * np.pi * 0.0026 for val in merged_resistance], marker=".", linestyle='-', c="black")

# Set the x-axis and y-axis labels for ax1
ax3.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax3.set_ylabel('Resistivity (µΩ.m)', size=20, labelpad=10)
ax3.set_title("Heat Up", size=16)
# margin: left,bottom,right,top
ax3.set_position([0.15, .17, 0.84, 0.77])
# Set the font size for x-ticks and y-ticks for ax1
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(72,150)
ax3.set_ylim(-3,25)
ax3.grid(True, linestyle='--', alpha=0.7)
ax3.figure.savefig("resistivity_up_with_line.svg", format = "svg")


# Create the inset plot of resistivity
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot the whole curve in the main plot
ax1.plot(merged_temperature, [val * 1000 for val in merged_resistance], marker=".",linestyle="-", c="black")
ax1.set_xlabel('Temperature (K)', size=20, labelpad=10)
ax1.set_ylabel('Resistance (mΩ)', size=20, labelpad=10)
ax1.set_title("Heat Up", size=16)
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
plt.savefig("resistance_up_inset.svg", format="svg")