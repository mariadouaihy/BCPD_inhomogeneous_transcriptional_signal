# BCPD_inhomogeneous_transcriptional_signal
Code associated with the article "Dissecting the dynamics of coordinated active transcriptional repression in a multicellular organism"


The BCPD method introduced in  Adams & MacKay 2007: is used to determine the onset of repression of a transcriptional signal imaged using the MS2-MCP system. The repression time for each singal independently is determed by the probability of having a change point denoting a sudden change in the parameters that generate the data

The code is organised as following: 
Input files: excel sheets with each column is a transcriptional intensity from one nuclei. The movies for the same sheets can be stored in different excel sheets in the same folder. The excel sheets should have the same formal as the example given in ./Data/
Output folder: ./Results/images/ contains the images of the nuclei traces per movies
               ./Results/npzFile/ contains the .npz format of the excel sheets
               ./Results/CP_0.6xMax_K_2/ contains the following with respect to each percentage, smoothing kernel given:
                       ./Results/CP_0.6xMax_K_2/full_signal/ contains the whole signal with the change point indicated by a change of color in the signal
                       ./Results/CP_0.6xMax_K_2/activation/ contains only the active part of the signal
                       ./Results/CP_0.6xMax_K_2/repressed/  contains only the repressed part of the signal
                       ./Results/CP_0.6xMax_K_2/T0_CP_0.6_K_2.xlsx contains the time point where the change point was found for each nuclei
                       ./Results/CP_0.6xMax_K_2/T0_cdf_1.pdf is the cumulative distribution function of the time into repression from the different nuclei.

The parameters needed are the following:
inputpath: path to the excel sheets
data_type: if we wan t to add a specific name to the output folder otherwise it's generate automatically
extension: extension type of the data (either xlsx or xls)
outputpath: output path where the Results folder will be created  
percentage: time point after percentage of max intensity
kernel:  kernel of smoothing the data
filtered: 1 if we want to remove nuclei that didn't go properly into repression from the statistics
p: probability of having a changepoint >=p

""" parameters of the signal and the movie """

FrameLen: frame length in seconds or time resolution
retention: in seconds
Polym_speed: Polymerase speed'
TaillePreMarq: length of mRNA until we reach the beginning of the MS2 sequence in bp
TailleSeqMarq: length of MS2 sequence in bp
TaillePostMarq: length of mRNA from last loop of MS2 until the polymerase leaves the TS
EspaceInterPolyMin: minimal distance between polymerase (in bp)
Intensity_for_1_Polym: calibration factor ( 1 is the data is already divided by the calibration)  
