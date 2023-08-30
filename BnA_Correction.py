### this py files implements the functions from Wagih et al R package in python
# see https://github.com/boonelab/sgatools
# Baryshnikova et al doi:10.1038/nmeth.1534
# Wagih et al doi:10.1093/nar/gkt400

import numpy as np
import pandas as pd
from tqdm import tqdm
from loess.loess_1d import  loess_1d
from matplotlib import pyplot as plt



pd.options.mode.chained_assignment = None # gets rid of the chain assignment warnings


##### N1 normalisation function ################
def N1_normalisation(DF, metrics ,default_cols=['gene','plate','row','column'], p_format=[32,48], default_overall_median= None):
    '''# Plate normalization: brings all plates to one same scale
    default_overall_median: dict storing {metric: screen median}, if None is calculated given the DF'''
    
    #check functions does job
    check=0
    
    #used to get the median of the 60% center colonies
    upper = 0.8
    lower = 0.2
    #we create an empty DF in which we will store the corrected metrics
    cols= default_cols
    [cols.append(m ) for m in metrics]
    DF_corrected = DF[cols]
    for metric in metrics: #['mean_intensity.24.YUGN','area.24.YUGN']:
        #if screen median is not provided we compute it here:
        if default_overall_median == None:
            i_select=  DF_corrected[metric].isna() == 0
            screen_median = DF_corrected[i_select][metric].median()
        else:
            
            screen_median = default_overall_median[metric]

    #create empty column:
        DF_corrected['corrected_' + metric]=  np.nan
       #> print(metric)
        
        for p in tqdm( np.unique(DF_corrected.plate)):
            # select the plate:
            #print(p,(DF_corrected[metric].isna() ==0).sum(),(DF_corrected[metric].isna() ==0).sum())
            i_select= (DF_corrected.plate == p) & (DF_corrected[metric].isna() ==0)
            
            #redefine upper and lower bounds, computes center 60% median, computes plate scaling factor
            upper = DF_corrected[i_select][metric].quantile(q=0.8)
            lower = DF_corrected[i_select][metric].quantile(q=0.2)
            i_select_60 = (DF_corrected[i_select][metric] >=lower ) & (DF_corrected[i_select][metric] <= upper)
            
            plate_median =DF_corrected[i_select][metric][i_select_60].median()
            scaling_factor = screen_median/plate_median
            #normalize plate in DF:
            DF_corrected.loc[i_select,'corrected_' + metric]=  DF_corrected.loc[i_select,metric]*scaling_factor

            # check
            if check:
                plt.title('plate ' +str(p)+' orange= scaled')
                plt.hist(DF_corrected[i_select][metric])
                plt.axvline(x =plate_median)
                plt.hist(DF_corrected[i_select]['corrected_' + metric], alpha=0.5)
                plt.axvline(x =screen_median, color='orange')
                plt.xlim([screen_median-1*screen_median, screen_median+1*screen_median])
                plt.show()

    return DF_corrected

##### N2 normalisation function ################
## Helper functions

def N2_spatial_normalisation_1plate(mat, padding_type = 'zeros' ):

    check=0 # for debuggin
    #Fill NA with a placeholder (mean of all colonies on plate) 
    i_nan = np.isnan(mat)
    
    if check :
        #define image boudaries for plotting on same scale
        vmin=  np.mean(mat[i_nan ==0])-2*np.std(mat[i_nan ==0])
        vmax= np.mean(mat[i_nan ==0])+2*np.std(mat[i_nan ==0])
        plt.imshow(mat, vmin=vmin,vmax=vmax)
        plt.show()
    mat[i_nan]= mat[i_nan==0].mean()
    #define Gaussian Kernel
    kernel_Gauss = fgaussian(7,2)
    #define average filter Kernel
    kernel_average = faverage(9)
  # Fill in NA with smoothed version of neighbors using gaussian blur
    mat_Gauss =  apply_filter(mat, kernel_Gauss , padding_type = 'zeros')
    mat[i_nan] = mat_Gauss[i_nan]

  # Apply median filter then average filter
    filtered = medianfilter2d(mat, 7, padding_type='replicate')
    filtered = apply_filter(filtered, kernel_average, 'replicate')
    # Subtract the mean of the filtered data from the filtered data
    filtered = filtered / np.mean(filtered)
  # Subtract filtered - mean from 
    normalized_mat = mat / filtered
    # put nan back into positions:
    normalized_mat[i_nan] = np.nan
    
    if check :
        plt.imshow(normalized_mat, vmin=vmin,vmax=vmax)
        plt.show()
        
    return normalized_mat

def medianfilter2d(mat, dim, padding_type = 'zeros'):
    '''mat: plate to filter
    dim: kernel size
    padding type can be 'zeros' anything else is edge (see np.pad) '''
    
    assert dim%2==1, "kernel dimension must be odd"
    check= 0
    r = np.zeros(mat.shape)
    
    margin_size = int((dim-1)/2)

    if padding_type == 'zeros':
        padded_mat = np.pad(mat, margin_size, mode='constant' )
    else:
        padded_mat = np.pad(mat, margin_size, mode='edge')
       
    filtered_mat = np.zeros(padded_mat.shape).astype(float)
    for r in range(margin_size, padded_mat.shape[0]-margin_size):
        for c in range(margin_size, padded_mat.shape[1]-margin_size):
  
            kernel = padded_mat[r-margin_size:r+1+margin_size,c-margin_size:c+1+margin_size]

            kernel = kernel[np.isnan(kernel)==0]
            filtered_mat[r,c]= np.median(kernel)
    if check:
        plt.imshow(padded_mat)
        plt.show()
        plt.imshow(filtered_mat[margin_size:-margin_size,margin_size:-margin_size])
        plt.show()
    
    return filtered_mat[margin_size:-margin_size,margin_size:-margin_size]

def fgaussian(x,sigma):
    check= 0
    x = np.arange(1,x+1)
    mat = np.empty((len(x),len(x)))
    for r in range(len(x)):
        for c in range(len(x)):
            n1=x[r]
            n2=x[c]
             # I think I spotted a mistake in Wagih's code
            #they define the kernel with -((n1)**2+(n)**2) when it should be: -((n1-np.mean(x))**2+(n2-np.mean(x))**2)
            mat[r,c] = np.exp(-((n1-np.mean(x))**2+(n2-np.mean(x))**2)/(2*sigma**2))
    if check:
        plt.imshow(mat/mat.sum())
    return mat/mat.sum()

def faverage(size):
    x = 1/(size*size)
    return np.repeat(x,size*size).reshape(size,size)


def apply_filter(mat, kernel, padding_type = 'zeros'):
    '''padding type can be 'zeros' anything else is edge (see np.pad) '''
    
    assert kernel.shape[0]%2==1 &kernel.shape[1]%2==1, "kernel dimension must be odd"
    assert kernel.shape[0] == kernel.shape[1], "kernel must be a square"
    
    check= 0 #for debuggin 

    margin_size = int((kernel.shape[0]-1)/2)

    if padding_type == 'zeros':
        padded_mat = np.pad(mat, margin_size, mode='constant' )
    else:
        padded_mat = np.pad(mat, margin_size, mode='edge')
       
    filtered_mat = np.zeros(padded_mat.shape).astype(float)
    for r in range(margin_size, padded_mat.shape[0]-margin_size):
        for c in range(margin_size, padded_mat.shape[1]-margin_size):
            mat_window = padded_mat[r-margin_size:r+1+margin_size,c-margin_size:c+1+margin_size]
            filtered_mat[r,c]= np.sum(mat_window * kernel)
    if check:
        plt.imshow(padded_mat)
        plt.show()
        plt.imshow(filtered_mat)
        plt.show()
    #return unpadded filetered array:
    filtered_mat=filtered_mat[margin_size: -margin_size,margin_size:-margin_size]
    return( filtered_mat)



def N2_spatial_normalisation(DF, metrics ,default_cols=['gene','plate','row','column'], p_format=[32,48]):
    '''# Spatial normalization: normalizes any gradient effect on the plate via median smoothing'''
     #check functions does job
    check=0

    #we create an empty DF in which we will store the corrected metrics
    cols= default_cols
    [cols.append(m ) for m in metrics]
    DF_corrected = DF[cols]
    for metric in metrics: #['mean_intensity.24.YUGN','area.24.YUGN']:
    #create empty column:
        DF_corrected[ metric]=  np.nan
        for p in tqdm( np.unique(DF_corrected.plate)):
            # select the plate:
            i_select= (DF_corrected.plate == p) 
            
            plate = DF[i_select][metric].values.reshape(p_format)
            normalised_plate = N2_spatial_normalisation_1plate(plate.copy(),padding_type='replicate')
            # we reshape the plate into its original shape
            DF_corrected.loc[i_select, metric]= normalised_plate.transpose(0,1).flatten()
            
            if check:
                #check that the reformatting of the plate matches that of DF
                n_equals= (plate.transpose(0,1).flatten() == DF[i_select][metric].values).sum()
                n_nan = np.isnan( DF[i_select][metric].values).sum()
                print(n_equals+n_nan)
                
                fig, axs= plt.subplots(1,2)
                axs = axs.flatten()
                i_nan = np.isnan(plate)
                vmin=  np.mean(plate[i_nan ==0])-2*np.std(plate[i_nan ==0])
                vmax= np.mean(plate[i_nan ==0])+2*np.std(plate[i_nan ==0])
                axs[0].imshow(plate)
                axs[1].imshow(normalised_plate)
                plt.show()
    return DF_corrected
    

    
##### (F3) Lowess Smoothing ####################
    
#Baryshnikova 2010:'trends across or down the plate should be relatively smooth'
#PC: it seems that they do not apply the LOWESS filter to each column. Rather they sort all colonies on the plate 
#according to their row/ column:

#for row smoothing
#x= [1,1...1,2,2,...,2,3,3,...,3...]
#y = [r1c1, r1,c2,...,r32,c1,r2c1,r2c2,...,r32c2,,r3c1,r3c2,...,r32c2..]
#at each row position, this gives a cloud of ncol points. they then apply the filter with a span of 6 rows

def N3_rowcol_Lowess(  DF,metrics,default_cols=['gene','plate','row','column'], p_format=[32,48]):
    '''we multiply each colony value by r_lowess_mean/r_lowess and by c_lowess_mean/c_lowess'''
    check=0
    cols= default_cols
    [cols.append(m ) for m in metrics]
    DF_corrected = DF[cols]
    for metric in metrics: #['mean_intensity.24.YUGN','area.24.YUGN']:
        #create empty column:
        DF_corrected[ metric]=  np.nan

        for p in tqdm( np.unique(DF_corrected.plate)):
            # select the plate:
            i_select= (DF.plate == p) & (DF[metric].isna()==0) 

            ## apply row wise LOWESS
            x = DF[i_select]['row'].values
            y = DF[i_select][metric].values
            if check:
                y_before = DF[DF.plate == p][metric].values
            i_sort = np.argsort(x)
            _, r_lowess, _ = loess_1d(x, y, frac=0.09)
            r_lowess_mean = np.mean(r_lowess)
            y_row_smoothed =  (y/r_lowess )*r_lowess_mean
            DF.at[i_select,metric] = y_row_smoothed

            ## apply column wise LOWESS
            x = DF[i_select]['column'].values
            y = DF[i_select][metric].values
            i_sort = np.argsort(x)
            _, c_lowess, _ = loess_1d(x, y, frac=0.09)
            c_lowess_mean = np.mean(c_lowess)
            y_rowcol_smoothed =  (y/c_lowess )*c_lowess_mean
            DF.at[i_select,metric] = y_rowcol_smoothed
            if check:
                
                fig, axs= plt.subplots(1,2)
                axs = axs.flatten()
                i_nan = np.isnan(y_before)
                vmin=  np.mean(y_before[i_nan ==0])-2*np.std(y_before[i_nan ==0])
                vmax= np.mean(y_before[i_nan ==0])+2*np.std(y_before[i_nan ==0])
                axs[0].set_title(metric)
                axs[0].imshow(y_before.reshape(32,48), vmin=vmin, vmax=vmax)
                i_select= (DF.plate == p) 
                axs[1].imshow(DF[i_select][metric].values.reshape(32,48), vmin=vmin, vmax=vmax)
                plt.show()

    return DF               

##### (F3) Jackknife filter ####################

def F3_jacknife_filter(DF,metrics,default_cols=['gene','plate','row','column'], p_format=[32,48]):
    '''applies JK filter to get rid of outliers
    JK is a leav one out method. We compute the std of quadruplicate colonies leaving 1 out: std_i.
    if std_i is very different from the mean of the 4 std_i (std_4) we get rid of the point that was left out.
    very differnt is define as greater than 0.9std_4'''

    check=0
    cols= default_cols
    [cols.append(m ) for m in metrics]
    DF_corrected = DF[cols]

    
    for metric in metrics: #['mean_intensity.24.YUGN','area.24.YUGN']:
        #create empty column:
        # select the plate:
        for p in tqdm( np.unique(DF_corrected.plate)):
            i_select= (DF_corrected.plate == p) & (DF_corrected[metric].isna()==0) 
            # in a few instances
            genes = np.unique(DF_corrected[i_select].gene)
            #we run through the genes in the plate, we assume that quadruplicates are placed side by side
            # we assume that empty positions for the metric have been NaNed
            for r in range(int(p_format[0]/2)):
                for c in range(int(p_format[1]/2)):
                    i_gene = i_select & (DF_corrected.row.isin([r*2,r*2+1])) & ( DF_corrected.column.isin([c*2,c*2+1]))
                    genes = DF_corrected[i_gene].gene.values
                    assert len(np.unique(genes))<2, ' quadruplicates of different genes'
                    #if any of the below conditions is not met we wont kife gene
                    knife_gene =  len(genes)>2 # we have 2 or less colonies
                    knife_gene = knife_gene & (False if len(genes)==0 else np.unique(genes) != '0') #the position shouldnt be empty
                    if knife_gene :
                        y = DF_corrected[i_gene][metric].values
                        #leave one out matrix
                        loo_matrix =  (np.identity(len(y))==0).astype(int) 

                        loo_matrix = loo_matrix*np.repeat(y.reshape(1,-1),len(y), axis=0) 
                        
                        #now the loo_matrix has replaced each looed value with 0. This will be taken to compute the std.
                        #we remove these 0s in the diagonal of loo_matrix with i_skip:
                        i_skip = [np.arange(len(y))!=i for i in np.arange(len(y))]
                        i_skip = np.vstack(i_skip)
                        loo_matrix = loo_matrix[i_skip].reshape(len(y),-1)
                        
                        #leave one out STDs
                        loo_std = np.std(loo_matrix , axis=1)
                        mean_loo_std = loo_std.mean()
                        # find outlier indexes (apply mean_loo_std – loo_std > mean_loo_std*(2/N) as in Baryshnikova 2010)
                        i_goodbye = np.abs(mean_loo_std - loo_std) > mean_loo_std*2/len(y)
 
                        #Wagih remove colonies that contribute more than 90% of total variance
                       #loo_var = np.var(loo_matrix , axis=1)
                        #total_var = np.var(y)
            
    

                        if check ==1:
                            if i_goodbye.sum()>0:
                                print(np.unique(genes) )
                                print('before knife ',y)
                                print(DF_corrected[i_gene][metric])
                        #replace the knifed outliers in y
                        y[i_goodbye] = np.nan
                        #replace the values in the main DF:
                        DF_corrected.at[i_gene, metric] = y

                        if check ==1:
                            if i_goodbye.sum()>0:
                                print('after knife ',y)
                                print(DF_corrected[i_gene][metric])
    return DF_corrected

