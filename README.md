# HLW_Digital_Twin
Time-resolved Digital Twin of Novel Cementitious Buffer Materials in High-Level Nuclear Waste Package

In order to make full use of the code, CrunchTope-v2.0.1 (the most recent release of the Crunchflow software) must also be installed. It is available at https://github.com/CISteefel/CrunchTope/releases/tag/v2.0.1.
Before running “Digital_twin.py,” update the simulation conditions using the CSV files in the “Input” directory.

The amount of each tracked isotope at time zero is entered into the “Initial_Isotope_inventory.csv” using the units of grams per canister.

The time of flight for the first reflection from two separate measurements of canister thickness (in units of seconds), canister wall temperature (in °C), and time between measurements (in years), as well as the time since time zero for the more recent of the two ultrasonic measurements (in years) are updated in “UT_data.csv.” 

Additional scenario assumptions are updated in “groundwater_and_dissolution.csv.” These include the minimum canister thickness needed to maintain canister integrity (in m), the flow rate of groundwater through a ruptured canister (in L/min), the assumed cross-section of the canister rupture (in m2), and the fractional leach rate for Am (fraction of Am in canister released per day).

Once the above files have been updated, run the “Digital_twin.py” script. This will produce a report with the predicted time to canister failure and expected isotope inventory at the time of rupture, which is saved to the “Output” directory along with plots of the isotope inventory and activity over time. It will also produce a Crunchflow input file. To complete the simulation, launch CrunchTope.exe from the “RT_directory.” This will simulate the uptake of Am on the AAM buffer. This will produce a report (time_series35.out) of the total Am concentration as a function of time since the canister ruptured at the AAM buffer repository boundary (35 cm from the canister).
