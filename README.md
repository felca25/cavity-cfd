# cavity-cfd

# This is a CFD solving for the classic lid-driven cavity problem
---

This code simulates the classic lid-driven cavity laminar flow problem from fluid dynamics.

The main control file is the "main.py file"

The solver is in the "cavity.py" file

The numerical functions used to calculate the parameters in the solver are located in the "functions.py" file

The file "save_data.py" has functions to stores data from any time_step choosen at the main function in .txt files in their respective folders, as it is more convenient and RAM saving to store data like so. An example of output can be seen at the "cavity-results" folder.

The "plot.py" file contains functions to plot the data from the .txt files.
