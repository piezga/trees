
import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.cm as mcolmap
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from itertools import chain
import scipy.sparse.csgraph as grp    
import scipy.io
import networkx as nx
import sys
import os
#sys.path.append(os.path.dirname(os.path.abspath("Crossvalidation")))
from sklearn.covariance import ShrunkCovariance
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.base import BaseEstimator
from sklearn.decomposition import PCA, FactorAnalysis
import pickle
import math
import estimators




def set_box_color(bp, color):
    '''
    sets the color of a boxplot
    '''
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)


# 'data' is TxN
def standardize_series(Xtrain,Xtest):
    avg = np.mean(Xtrain, axis=0)
    Xtrain -= avg 
    std = np.std(Xtrain, axis=0)
    Xtrain /= std
    
    Xtest -= avg
    Xtest /= std


def percolation_analysis(matrix,ndivisions=1000):
    M=np.copy(matrix)

    N,N=np.shape(M)
    M-=np.eye(N)*np.diag(M)
    M=np.abs(M)

    #indices = matrix < 0.
    #matrix[indices] = 0


    epsilon=1.0E-8
    
    maxv = np.max(M)
    delta = maxv/ndivisions

    n_ccs = np.zeros(ndivisions)
    number_cutted_links = np.zeros(ndivisions)

    
    for j in range(ndivisions):
        mask_indices = M < j*delta + epsilon
        M[mask_indices]=0.
        n_cc = grp.connected_components(M, directed=False, return_labels=False)
        n_ccs[j] = n_cc
        number_cutted_links[j] = np.sum(mask_indices)

    ncouples = N*(N-1.)
    fraction_surviving_links = 1. - (number_cutted_links-N)/ncouples
    return delta*np.arange(ndivisions) , n_ccs/N , fraction_surviving_links

def perc_analysis(matrix):
    matrix=np.abs(matrix)
    dic={}
    N,N=np.shape(matrix)
    dic[0]=0
    soglia_=0
            
    lista=np.unique(matrix) #ordina gli elementi della matrice, senza ripetizioni
    S=len(lista)-1
    i=S
    for thisno in range(1,N):
        M=np.copy(matrix)
        #lista=lista[S-i:]
        S=len(lista)-1

        if len(lista)>1000:  #parto da qui se la lista è più lunga di 1000
            for i in range(int(len(lista)/1000)): #scorre dal basso all'alto, ad intervalli di 1000

                j=i*1000
                ele= lista[j]

                indici=np.where( M < ele )
                M[indici]=0   
                cc=grp.connected_components(M, directed=False, return_labels=False)

                if  cc>thisno:

                    S=j
                    #print (i,cc, ele,'inverti') #indice 
                    break

        if len(lista)>100: #parto da qui se la lista è più lunga di 100
            for i in range(int(S/100)): #scorre dall'alto al basso, a passi di 100
                M=np.copy(matrix)
                ele= lista[S-i*100]
            #print (ele)
                indici=np.where( M < ele )
                M[indici]=0   
                cc=grp.connected_components(M, directed=False, return_labels=False)
                #print (cc)
                if  cc<thisno+1:
                    S=S-i*100
                    break

        if len(lista)>10: #parto da qui se la lista è più lunga di 10
            for i in range(int(S/10)): #scorre dal basso all'alto, a passi di 10
                M=np.copy(matrix)
                ele= lista[S+i*10]
                #print (ele)
                indici=np.where( M < ele )
                M[indici]=0   
                cc=grp.connected_components(M, directed=False, return_labels=False)
                #print (cc)
                if  cc>=thisno+1:

                    S=S+i*10
                    #print ('inverti')
                    break

        #questa la esegue sempre    
        for i in range(int(S)): #scorre dall'alto al basso, eleemnto per elemento 
            M=np.copy(matrix)
            ele= lista[S-i]
        #print ('ele',ele)
            indici=np.where( M < ele )
            M[indici]=0   
            cc=grp.connected_components(M, directed=False, return_labels=False)
            #print (cc)
            if  cc<=thisno:
                soglia_=ele
                
                #print ('fine',soglia_J)
                break
        dic[thisno]=soglia_
    return dic    



def regularize(meth,sub,a,**kwargs):
    if meth=='LS':
        return linear_shrinkage(np.corrcoef(sub.T),a)
    elif meth=='GS':
        return estimators.get_covariance(estimators.Shrinkage_biasedmatrix(alpha=a,M0=kwargs.get('T', None)).fit(sub)) 
    elif meth=='PCA':
        return estimators.get_covariance(PCA(svd_solver='full',n_components=a).fit(sub))
    elif meth=='consPCA':
        return  estimators.get_covariance(estimators.ConservativePCA(p=a).fit(sub))
    

def loadfile(file_,folder=''): #l'argomento deve essere una stringa (metti le virgolette)

    filename=folder+'file_'+file_+'.pkl'
    file = open(filename, "rb")
    loaded=pickle.load(file)
    file.close()
    return loaded



def mediafisher(M,thisax=0): 
    '''
    Inputs: M, to be averaged along axis thisas with the Fisher transform method 
    Output: Fisher averaged M
    '''
    N=np.copy(M)
    for C in N:
        C-=np.eye(len(C))
    return np.tanh(np.average(np.arctanh(N),axis=thisax))+ np.eye(len(N[0]))

def linear_shrinkage(M,alfa):
    '''
    linear shrinkage of squared matrix M with param alpha (0= no shrinkage, 1= identity matrix)
    '''
    return (1-alfa)*M+alfa*np.eye(len(M))

def C_media_shr_fisher(data,a): 
    '''
    Inputs: data-> subjects data to be regularized and averaged
            a-> regularization parameter(s)
    
    The regularized subjects are averaged with the Fisher transform 
    '''
    a = np.asarray(a)
    if a.ndim == 0:
        a = a[None]  # Makes x 1D
        scalar_input = True
        #print (np.squeeze(a))
        return(mediafisher([linear_shrinkage(np.corrcoef(sub.T),a)  for sub in data]))  
    else:
        return(mediafisher([linear_shrinkage(np.corrcoef(sub.T),a_s)  for sub,a_s in zip(data,a)]))

def standarise(M,thisax=0): 
    '''
    standardizes M along given axis
    '''
    Ms=np.copy(M)
    Ms-=np.average(Ms,axis=thisax)
    Ms/=np.std(Ms,axis=thisax)
    return Ms    

def standarise_matrix(M): 
    Mc=np.copy(M)
    Mc/=np.sqrt(np.outer(np.diag(Mc),np.diag(Mc)))
    return Mc    




def percolazione(matrix,graph=True,abs=True): #restituisce il grafo della matrice M sogliato secondo tale criterio
    M=np.copy(matrix)

    s=soglia_perc(np.abs(M))
    ind=np.where(np.abs(M)<s)

    M[ind]=0
    if graph:
        return  nx.from_numpy_matrix(np.abs(M))
    elif abs:
        return np.abs(M)
    else:
        return M

def soglia_perc(matrix):
    soglia_=0
    M=np.abs(np.copy(matrix))
    lista=np.unique(M) #ordina gli elementi della matrice, senza ripetizioni
    S=len(lista)-1
    
    if len(lista)>1000:  #parto da qui se la lista è più lunga di 1000
        for i in range(int(len(lista)/1000)): #scorre dal basso all'alto, ad intervalli di 1000

            j=i*1000
            ele= lista[j]

            indici=np.where( M < ele )
            M[indici]=0   
            cc=grp.connected_components(M, directed=False, return_labels=False)

            if  cc>1:

                S=j
                #print (i,cc, ele,'inverti') #indice 
                break
    
    if len(lista)>100: #parto da qui se la lista è più lunga di 100
        for i in range(int(S/100)): #scorre dall'alto al basso, a passi di 100
            M=np.abs(np.copy(matrix))
            ele= lista[S-i*100]
        #print (ele)
            indici=np.where( M < ele )
            M[indici]=0   
            cc=grp.connected_components(M, directed=False, return_labels=False)
            #print (cc)
            if  cc<2:
                S=S-i*100
                break
    
    if len(lista)>10: #parto da qui se la lista è più lunga di 10
        for i in range(int(S/10)): #scorre dal basso all'alto, a passi di 10
            M=np.abs(np.copy(matrix))
            ele= lista[S+i*10]
            #print (ele)
            indici=np.where( M < ele )
            M[indici]=0   
            cc=grp.connected_components(M, directed=False, return_labels=False)
            #print (cc)
            if  cc>=2:

                S=S+i*10
                #print ('inverti')
                break
    
    #questa la esegue sempre    
    for i in range(int(S)): #scorre dall'alto al basso, eleemnto per elemento 
        M=np.abs(np.copy(matrix))
        ele= lista[S-i]
    #print ('ele',ele)
        indici=np.where( M < ele )
        M[indici]=0   
        cc=grp.connected_components(M, directed=False, return_labels=False)
        #print (cc)
        if  cc<=1:
            soglia_=ele
            #print ('fine',soglia_J)
            break
    return soglia_


def stdMatrix(C):
    M=np.copy(C)
    aux=np.sqrt(np.outer(np.diag(M),np.diag(M)))
    return M/aux

def plot2matrices(mat1,mat2,filename=None,removediagonal=True,xlabel=None,ylabel=None,title=None,sci=True,logscale=False):

    if logscale:
        mat1_copia = (np.log10(mat1)).copy()
        mat2_copia = (np.log10(mat2)).copy()
    else:
        mat1_copia = mat1.copy()
        mat2_copia = mat2.copy()
    D=len(mat1)
    mat3=np.ones((D,D)) #matrice di tutti uni di taglia D
    for i in range(D):
        for j in range(D):
            if i<j:
                mat3[i,j]=mat1_copia[i,j] 
            else:
                mat3[i,j]=mat2_copia[i,j] #scrive una matrice con triangolo superiore uguale a mat1 e triangolo inferiore + diagonale uguale a mat2
    if removediagonal:
        mat3-=np.eye(D)*np.diag(mat3) 
    plt.matshow(mat3) 
    
    if sci:
        plt.colorbar(format='%.0e')
    else:
        plt.colorbar()
        
    #plt.xticks(range(D))
    #plt.yticks(range(D))
    if title is not None:
        plt.title(title,fontsize=28)
        plt.subplots_adjust(top=0.80)
        
    if xlabel is not None:
        plt.xlabel(xlabel,fontsize=30)
    if ylabel is not None:
        plt.ylabel(ylabel,fontsize=30)
    if filename is not None:
        plt.savefig(filename)
    plt.show()

def spec_ent_fast(G, betas,normalize=True):
    G.remove_edges_from(nx.selfloop_edges(G))
    spect=nx.laplacian_spectrum(G)
    if normalize:
        spect*=len(spect)/np.sum(spect)
        
    expltau=np.exp(np.outer(-1*betas,spect)) #e^-tau lambda_k (T righe, N colonne)
    Z=np.sum(expltau,axis=1) #vettore N
    return np.sum(expltau*((np.log(Z)+np.outer(betas,spect).T).T),axis=1)/Z# T righe, N colonne



def betasfunct(G,n=30,weight='weight'): 
    
    spect=np.array(sorted(nx.laplacian_spectrum(G,weight=weight)))
    
    ref= np.log(nx.number_of_nodes(G))

    #first gruess
    maxbeta=10/spect[1]
    minbeta=n/np.sum(spect)
    thislist=np.logspace(math.log10(minbeta),math.log10(maxbeta),num=n)
        
    #print (spect,taumin,taumax)

    
    #check on min  (not too small)
    a=(sum((thislist[3]*spect+np.log(sum(np.exp(-thislist[3]*spect))))*np.exp(-thislist[3]*spect)))/sum(np.exp(-thislist[3]*spect)) 
    while np.abs(a-ref)<10E-2: 
        minbeta+=((thislist[1]-thislist[0])) # aumento il minimo
        thislist=np.logspace(math.log10(minbeta),math.log10(maxbeta),num=n)
        a=(sum((thislist[3]*spect+np.log(sum(np.exp(-thislist[3]*spect))))*np.exp(-thislist[3]*spect)))/sum(np.exp(-thislist[3]*spect)) 

    
    #check on min  (not too big)       
    a=(sum((minbeta*spect+np.log(sum(np.exp(-minbeta*spect))))*np.exp(-minbeta*spect)))/sum(np.exp(-minbeta*spect)) 
        
    while np.abs(a-ref)>10E-3: 
        #print (minbeta,np.abs(a-ref))
        i=1/10
        minbeta-=((thislist[1]-thislist[0])/(i*10))
        #print (minbeta)
        if minbeta<=0:
            minbeta+=(thislist[1]-thislist[0])
            i=1
        thislist=np.logspace(math.log10(minbeta),math.log10(maxbeta),num=n)
        a=(sum((minbeta*spect+np.log(sum(np.exp(-minbeta*spect))))*np.exp(-minbeta*spect)))/sum(np.exp(-minbeta*spect)) 
        
     #check on max (not too small)        
    a=(sum((maxbeta*spect+np.log(sum(np.exp(-maxbeta*spect))))*np.exp(-maxbeta*spect)))/sum(np.exp(-maxbeta*spect)) 

    while a>10E-5:
       # print (a)
        maxbeta+=(thislist[-1]-thislist[-2])
        a=(sum((maxbeta*spect+np.log(sum(np.exp(-maxbeta*spect))))*np.exp(-maxbeta*spect)))/sum(np.exp(-maxbeta*spect)) 
        thislist=np.logspace(math.log10(minbeta),math.log10(maxbeta),num=n)
        
    
    return thislist

def binarizza(N,soglia):
    M=np.copy(N)
    indici=np.where(np.abs(M)<soglia)
    for el in zip(indici[0],indici[1]):    
        M[el]=0. 
    indici=np.where(np.abs(M)>0)
    for el in zip(indici[0],indici[1]):    
        M[el]=1. 
    
    return M

def Mdist(M1,M2):
    '''calcola la correlazione fra gli elementi delle matrici 1 e 2'''
    return np.sum(M1*M2)/np.sqrt(np.sum(M1**2)*np.sum(M2**2))
def Mdist_norm(M1,M2):
    '''
    calcola la correlazione fra gli elementi delle matrici 1 e 2, al netto dei link comuni
    '''
    return np.sum(M1*M2)/np.sqrt(np.sum((M1*binarizza(M2,0))**2)*np.sum((M2*binarizza(M1,0))**2))

