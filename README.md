# Analysis_for_Time_Resolution

This is time resolution of radio particle detector with ROOT. 
To do analysis, please install ROOT from https://root.cern.

This is analysis for time resolution.
The time resolution is the fluctuation of the response time of the detector, which is the time that the events within that time cannot be separated.
This analysis is for qdc and tdc data from particle and nuclear experiments.

In ana112.C Using the make class of ROOT's tree, obtain the values of pedestal and MIP, make corrections related to energy (slewing correction), and obtain the time resolution.
Be careful with this order.
Other analysis macros can be used to compare the time resolution obtained with ana112.C or to verify the conditions of the cut in ana112.C.

However, the index of the array of data to be read depends on my experimental environment, so if you use it, please modify it to fit your environment.
Please use run10000000027.root as sample data.
