import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import shutil
import zipfile
from scipy.io import loadmat

#from joblib import Parallel, delayed

### imput: index of the file
class readDataInFile:
    def __init__(self,DataFilePath,outputFolder,data_type, fParam, extension='.xlsx'):
        dataFile='npzFile'+data_type
        PyFilePath = os.path.join(outputFolder,dataFile)
        ### where the images will be stored
        DataFilePath0 = os.path.join(outputFolder,'images')
        ### creating the folder
        if not os.path.exists(PyFilePath):
            os.mkdir(PyFilePath) 
 
        if not os.path.exists(DataFilePath0):
            os.mkdir(DataFilePath0) 
            
        self.PyFilePath = PyFilePath
        self.DataFilePath0=DataFilePath0
        self.DataFilePath=DataFilePath
        
        self.lists, self.nexp = self.ListingFiles(self.DataFilePath, extension=extension)    
        self.file_name_list=self.lists
        for data_i in range(self.nexp):
            if '.~' not in self.file_name_list[data_i]:
                self.a=self.read_data(data_i, fParam, outputFolder, extension)
    def read_data(self,data_i, fParam, outputFolder,extension):

        if '.npz' in fParam:
            deconParameters=np.load(fParam)
    
        elif '.mat' in fParam:
            deconParameters=loadmat(fParam)
    
        Polym_speed = deconParameters['Polym_speed'] 
        TailleSeqMarq = deconParameters['TailleSeqMarq']
        TaillePostMarq = deconParameters['TaillePostMarq']
        EspaceInterPolyMin = deconParameters['EspaceInterPolyMin']
        FrameLen = deconParameters['FrameLen']
        Intensity_for_1_Polym = deconParameters['Intensity_for_1_Polym']
        FreqEchImg = deconParameters['FreqEchImg']
        DureeSignal = deconParameters['DureeSignal']
    
        file_name_list=self.file_name_list

        DataFileName = self.DataFilePath + file_name_list[data_i] #path to the specified file
        print('File '+ str(data_i+1) + ': ' + file_name_list[data_i])
        if not zipfile.is_zipfile(DataFileName):
            DataExp =pd.read_excel(DataFileName,usecols=np.arange(4,300), skiprows = range(0,1)).to_numpy() #, engine='openpyxl').to_numpy()
            Spot_names=pd.read_excel(DataFileName,usecols=np.arange(4,len(DataExp[0])+4),skipfooter =len(DataExp)+1).columns.tolist() 
        else:
            DataExp =pd.read_excel(DataFileName,usecols=np.arange(4,300), skiprows = range(0,1), engine='openpyxl').to_numpy() 
            Spot_names=pd.read_excel(DataFileName,usecols=np.arange(4,len(DataExp[0])+4),skipfooter =len(DataExp)+1, engine='openpyxl').columns.tolist()

        
        ## creating a folder to store the images
        dirname=file_name_list[data_i].replace( extension,'')
        
        
        
        
        WriteTo=self.DataFilePath0+'/'+dirname+'_image'
 
        if not os.path.exists(WriteTo):
            os.mkdir(WriteTo) 
            
        ## file name and size of dataexp
        n=DataExp.shape
        fn='data_'+self.file_name_list[data_i].replace('_','').replace(extension,'')

        ## ploting the intensity spot for each nuclei in a 6*6 plot
        for iii in range(n[1]):

            ifig= math.floor((iii)/36)+1 
            h=plt.figure(ifig+data_i*10, figsize=(11.69,8.27))
            plt.subplot(6,6, (iii%36+1))
            plt.plot(np.arange(1,len(DataExp)+1), DataExp[:,iii], color='black', linewidth=0.1)
            h.subplots_adjust(hspace = 1, wspace=0.4)

            plt.title(Spot_names[iii].replace('_','-'))

            if (iii+1)%36==0 or (iii+1)== n[1]:

                figfile=WriteTo+'/fig'+str(ifig)+'.pdf'
                h.savefig(figfile)
                plt.close(h)
                
                
        ### removing zeros nuclei
        idx_zeros = np.where(np.sum(DataExp,axis=0)==0)[0] #avoid zero cells
        if sum(idx_zeros)!=0:
            idx_spots = np.delete(np.arange(DataExp.shape[1]), idx_zeros)
            removed_zeros_spot = np.array(Spot_names)[idx_zeros].tolist()
            Spot_names = np.array(Spot_names)[idx_spots].tolist()
        DataExp = np.delete(DataExp, idx_zeros, axis=1)
        
        
        ### removing nuclei that has an estimation of less then 1 polII per nc
        sd=DataExp.shape
        frame_num=sd[0] ### number of frames
        DureeSimu = frame_num*FrameLen  ### film duration in s
        DureeAnalysee = DureeSignal + DureeSimu ###(s)
    
        num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)) # maximal number of polymerase positions
        area = FreqEchImg/Polym_speed*Intensity_for_1_Polym*(TaillePostMarq+TailleSeqMarq/2)
        DataExpSmooth = np.minimum(np.round(np.sum(DataExp,axis = 0) / area), num_possible_poly)
        idx_zeros_smoothed = np.where(DataExpSmooth==0)[0]
        if sum(idx_zeros_smoothed)!=0:
            idx_spots = np.delete(np.arange(DataExp.shape[1]), idx_zeros_smoothed)
            removed_low_spot = np.array(Spot_names)[idx_zeros_smoothed].tolist()
            Spot_names = np.array(Spot_names)[idx_spots].tolist()
        DataExp = np.delete(DataExp, idx_zeros_smoothed, axis=1)
                
        ## storing removed nuclei
        f = open(outputFolder+'removed_nuclei', "a")
        f.write("\n "+ dirname +':')
        if sum(idx_zeros)!=0:
            f.write("\n \t  zeros nuclei: " + str(removed_zeros_spot))
        if sum(idx_zeros_smoothed)!=0:
            f.write("\n \t  nuclei with < 1 possible pol: " + str(removed_low_spot))
        f.close()
        
        #saving the data
        sd=DataExp.shape      
        Frames=np.arange(0,sd[0])
        Samples=np.arange(1,sd[1])
        fname=self.PyFilePath+fn+'.npz'
        self.DataExp=DataExp
        self.Frames=Frames
        self.Samples=Samples

        np.savez(fname,DataExp=DataExp,Samples=Samples,Frames=Frames)
        self.fname=fname
        return self.DataFilePath0

    def ListingFiles(self,DataFilePath,extension):
        self.DataFilePath=DataFilePath
        file_name_list = np.array(os.listdir(self.DataFilePath)) # list of subdirectories containing data from different genotypes
        nexp = len(file_name_list) # number of files inside he folder
        #  DataFilePath0='../images/' # place to store images
        #just to avoid erros in the opening of the files
        file_name_list_modified=np.array([])
        for i in range(nexp):
                if '._' in file_name_list[i] :
                    file_name_list[i]=file_name_list[i].replace('._','',1)
                    
        for i in range(nexp):
            if not '._' in file_name_list[i] and extension in file_name_list[i]:
                file_name_list_modified=np.append(file_name_list_modified,file_name_list[i])
                    

        file_name_list=file_name_list_modified            
        file_name_list=np.unique(file_name_list)
        nexp=len(file_name_list)
        return file_name_list, nexp
            
        
            
