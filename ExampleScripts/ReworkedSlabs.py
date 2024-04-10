from chrest.ablateData import AblateData
import matplotlib.pyplot as plt
import numpy as np
import math as math

#######################################################################################################################
# Specify the location where your hdf5 files are located:
File_Location = "/p/lustre2/ubchrest/post23Tests/test2D-dm0-pmma-2262275/domain/domain."
# Specify the number of hdf5 files in that folder:
NFiles = 246
# Specify the interval of files read reports:
Report_Number = 50
# Specify number of interpolation points
Interpolation_Number = 1000
# Define an x-coordinate (meters) where a slice should be taken:
x_slice = 0.07
# Choose an interpolation method:
method = "nearest"
# Select a timestep to view temperature data:
selected_timestep = 50
# Select a plot-type (XSlice, in-development)
Selected_Plot = "XSlice"

#######################################################################################################################
# Define functions
def initializing():
    # Initialization data
    print("Initializing")
    initialization_file = File_Location + "00000" + ".hdf5"
    initialization_data = AblateData(initialization_file)
    initialization_tempdata = initialization_data.get_field('aux_temperature')
    initialization_arrdata = initialization_tempdata[0]
    NCells = len(initialization_arrdata[0])
    NFileArray = np.arange(NFiles)
    Tempdata = np.zeros((NFiles, NCells))
    Presdata = np.zeros((NFiles, NCells))
    Timedata = np.zeros(NFiles)
    TempAllTimes = np.zeros((NFiles, Interpolation_Number))
    PresAllTimes = np.zeros((NFiles, Interpolation_Number))
    return NCells, NFileArray, Tempdata, Presdata, Timedata, TempAllTimes, PresAllTimes
def pulling_data():
    print("Pulling data from folder")
    for i in NFileArray:
        # Create file variable
        number = f'{(i):05}'
        file = File_Location + number + ".hdf5"
        if i % Report_Number == 0:
            # Report the file to user.
            print("Reading data from " + file)
        # Create the Ablate data object for each specific file.
        data = AblateData(file)
        # Read the temperature data for specific hdf5 file.
        Temperature = data.get_field('aux_temperature')
        # Read the pressure data for specific hdf5 file.
        Pressure = data.get_field('aux_pressure')
        # Read the time data.
        Times = data.times
        # Translate the read data into an array.
        Tempdata[i, :] = Temperature[0]
        Presdata[i, :] = Pressure[0]
        Timedata[i] = Times[0]
        # Create variables for cell data.
        Cells = data.cells
        CellCenterCoords = data.compute_cell_centers()
        YValues = CellCenterCoords[:, 1]
        XValues = CellCenterCoords[:, 0]
    print("All specified data read")
    return Cells, CellCenterCoords, Timedata, Presdata, Tempdata, YValues, XValues
def interpolation():
    print("Starting data interpolation")
    from scipy.interpolate import griddata as gd
    # Define the interpolation points.
    min_y = min(CellCenterCoords[:, 1])
    max_y = max(CellCenterCoords[:, 1])
    y = np.linspace(min_y, max_y, Interpolation_Number)
    x = np.ones(Interpolation_Number) * x_slice
    # Convert the points into Xi, the interpolation location coordinates.
    Xilist = [x, y]
    Xilist = np.transpose(Xilist)
    Xi = np.asarray(Xilist)
    # Interpolate the data onto selected Xi.
    for tidx in NFileArray:
        tidx = int(tidx)
        # Pressure data extraction
        Interpolated_Pres = gd(CellCenterCoords, Presdata[tidx, :], Xi, method=method)
        PresAllTimes[tidx] = Interpolated_Pres
        # Pressure data extraction
        Interpolated_Temp = gd(CellCenterCoords, Tempdata[tidx, :], Xi, method=method)
        T = Tempdata[tidx, :]
        TempAllTimes[tidx] = Interpolated_Temp
    print("Finished data interpolation")
    return x, y, TempAllTimes, Interpolated_Temp, PresAllTimes, Interpolated_Pres

#######################################################################################################################
# Pulling data from functions
if __name__ == "__main__":
    NCells, NFileArray, Tempdata, Presdata, Timedata, TempAllTimes, PresAllTimes = initializing()
    Cells, CellCenterCoords, Timedata, Presdata, Tempdata, YValues, XValues = pulling_data()
    x, y, TempAllTimes, Interpolated_Temp, PresAllTimes, Interpolated_Pres = interpolation()

#######################################################################################################################
# Plotting
    counter = np.linspace(0, Interpolation_Number-1, Interpolation_Number)

    if Selected_Plot == "Temp vs Time":
        plt.figure(1)
        fig, ax = plt.subplots(figsize=(10, 10))
        ColorControl = np.linspace(0, .75, len(counter))
        for i in counter:
            i = int(i)
            j = np.round(y[i], 3)
            ax.plot(Timedata, TempAllTimes[:, i], color=str(ColorControl[i]), lw=2) # , label="y = " + str(j)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temp (K)")
        #ax.legend(loc="upper right")
        plt.title("Temperature along sample points at x = .07 m")
        fig.show()

    if Selected_Plot == "XSlice":
        # Determine the location of the slab by seeing where pressure = 0
        CutoffPressure = np.transpose(np.asarray(PresAllTimes))
        index = []
        for i, counting in enumerate(CutoffPressure):
            if PresAllTimes[50, int(i)] == 0:
                index = np.append(index, i)
                if (i > 1) & (i - int(np.size(index)) > 1):
                    cut1 = int(np.size(index)) - 1
                    cut2 = i
                    break

        cropped_pressure = PresAllTimes[selected_timestep, cut1:cut2]
        cropped_temperature = TempAllTimes[selected_timestep, cut1:cut2]

        # Pressure vs X - Slice
        plt.figure(2)
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.plot(y * 100, PresAllTimes[selected_timestep, :] / 101325, color="k", lw=2)
        ax.set_xlabel("Y Coordinate (cm)", fontdict={'fontsize':16})
        ax.set_ylabel("Pressure (atm)", fontdict={'fontsize':16})
        plt.title("2D Pressure", fontdict={'fontsize':22}, fontweight='bold')
        plt.legend(["x = " + str(round(x_slice*100)) + "cm, t = " + str(round((Timedata[selected_timestep] * 1000), 3)) + "ms"],  prop={"size": 15})
        plt.xlim((y[cut1] * 100), (y[cut2] * 100)-.01)
        plt.ylim(np.min(cropped_pressure / 101325), round(np.max(np.asarray(cropped_pressure / 101325)), 2)+.005)
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        plt.savefig( '/p/lustre2/ubchrest/post23Tests/test2D-dm0-pmma-2262275/domain/IntPressure.png', bbox_inches='tight')
        fig.show()

        # Temperature vs X - Slice
        plt.figure(3)
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.plot(y*100, TempAllTimes[selected_timestep, :], color="k", lw=2)
        ax.set_xlabel("Y Coordinate (cm)", fontdict={'fontsize':16})
        ax.set_ylabel("Temperature (K)", fontdict={'fontsize':16})
        plt.title("2D Temperature", fontdict={'fontsize':22}, fontweight='bold')
        plt.legend(["x = " + str(round(x_slice*100)) + "cm, t = " + str(round((Timedata[selected_timestep] * 1000), 3)) + "ms"], prop={"size": 15})
        plt.xlim((y[cut1] * 100)+.01, (y[cut2] * 100))
        plt.ylim(math.ceil((np.min(cropped_temperature),2)[0]), math.ceil((np.max(np.asarray(cropped_temperature)), 2)[0])+100)
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        plt.savefig('/p/lustre2/ubchrest/post23Tests/test2D-dm0-pmma-2262275/domain/InterpolatedTemperature.png', bbox_inches='tight')
        fig.show()

        # Velocity Components vs X - Slice
        # Velocity Magnitude vs X - Slice

    print("Done")