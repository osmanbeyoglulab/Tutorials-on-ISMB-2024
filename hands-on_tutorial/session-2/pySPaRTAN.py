
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Xiaojun Ma
# Created Date: Mar 18 10:54:00 PDT 2020
# =============================================================================
"""
This script contains the major class of SPaRTAN model and its dependencies.

This script requires numpy, Scipy, matplotlib to be installed within the Python
environment you are running this script in

This script requires Cython modules present in the current directory

This file contains the following classes and functions

    class Timer: a class to convert time period in seconds to the format of h:m:s

    class pySPaRTAN: The major class for SPaRTAN, establishing an interaction matrix between
    surface proteins (P) and TFs (D) that predict target gene expression (Y).

    function normalize_column(): perform l2 normalization column-wize of given matrix

"""
import numpy as np
import pandas as pd
import cythKronPlus as krnP
import cythLeastR as leastR
import scipy.linalg
import time
import gc
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.model_selection import KFold

class Timer:
    """ a class to convert time in seconds to the format of h:m:s

    Methods
    -------
    def __init__(self):
        initiate a timer

    def restart(self):
        restart a timer

    def get_time_hhmmss(self):
        return the period = end_time - start_time in (h, m, s) format

    """

    def __init__(self):
        # initiate a timer
        self.start = time.time()

    def restart(self):
        # restart a timer
        self.start = time.time()

    def get_time_hhmmss(self):
        # return the period = end_time - start_time
        # in (h, m, s) format
        end = time.time()
        m, s = divmod(end - self.start, 60)
        h, m = divmod(m, 60)
        time_str = "%02d:%02d:%02d" % (h, m, s)
        return time_str


def normalize_column(A, T=0):
    """ perform l2 normalization column-wize of given matrix

    Parameters:
        A : the matrix that works on
        T : switch of column-wize and row-wize.
            T=0: column-wize
            T=1: row-wize

    """
    if (T == 0):
        return np.divide(A, np.sqrt(np.sum(A**2, 0)))
    else:
        At = np.transpose(A)
        return np.transpose(np.divide(At, np.sqrt(np.sum(At**2, 0))))


class pySPaRTAN:
    """
    The major class for SPaRTAN, establishing an interaction matrix between
    surface proteins (P) and TFs (D) that predicts target gene expression (Y).

    Methods
    -------
    fit(self, D, P, Y)
        cross validate and train a SPaRTAN model

    _fit(self, D, P, Y,lamda=None, rsL2=None,  spectrum=None )
         train a SPaRTAN model with specific hyperparameters
         
    _cv(D, P, Y)
         cross validation to determine the best hyperparameters
         
    _ar_model2w(self):
        converts a trained model to intermidiat vaiable W

    _ar_reconstruction(self, pred_test=None):
        reconstruction function

    predict(self, P_test=None):
        predict target gene expression

    get_corr(self, Y_pred, Y_test, plot=False):
        get the correlation between predicted Y_pred and Y_test

    get_W(self):
        get coefficient matrix

    get_projP(self, Y=None):
        get projected protein expression

    get_projTF(self, P=None):
        get projected TF activity

    """
    
    def __init__(self, lambdas=None, rsL2s=None,  spectrums=None, corrtype='pearson', n_fold=5):
        
       
        self.rsL2s = rsL2s
        self.lambdas = lambdas
        self.spectrums = spectrums
  
        self.corrtype = corrtype
        self.n_fold=n_fold
        
    def  _cv(self, D, P, Y):
        
      
       
        lensps = len(self.spectrums)
        lenlambdas = len(self.lambdas)
        lenrsL2s = len(self.rsL2s)
  
          
        corr_all = np.zeros((lensps, lenlambdas, lenrsL2s))
        
        # for a in range(0, lenspAs):
        for s in range(lensps) :
            for l in range(lenlambdas):
                for r in range(lenrsL2s):
                    print(
                        "cross validating  spectrum={}, lambda={}, rsL2={}".format(
                            self.spectrums[s], self.lambdas[l] ,self.rsL2s[r]
                        )
                    )
                    Y_pred_all = Y * 0


                    #added some parameters
                    kf = KFold(n_splits = self.n_fold, shuffle = True, random_state = 2)
                    result = next(kf.split(P), None)

                  
                    corrs = 0
                    for f in range(self.n_fold):

                         # split the data into train and test set
                        P_train, P_vali = P.iloc[result[0], :], P.iloc[result[1], :]
                        Y_train, Y_vali = Y.iloc[result[0], :], Y.iloc[result[1], :]

            
                        # train the model
                        self._fit(
                            D,
                            P_train,
                            Y_train,
                            lamda=self.lambdas[l],
                            rsL2=self.rsL2s[r],
                            spectrum=self.spectrums[s],
                        )

                        # get predicted value Y_pred  on P_vali
                        Y_pred = self.predict(P_vali)

                        # get the correlation bewteen Y_pred and Y_vali
                        Y_pred_all.iloc[result[1],:] = Y_pred
                        
                        corrs += self.get_pred_score(P_vali, Y_vali)


                    corr = corrs/self.n_fold
                    corr_all[s, l, r] = corr


        # outfile_corr = os.path.join(outpath,"cvPerform.txt")
        # output_corr(outfile_corr, corr_all)
        # logMem("Cross validation done", memfile)

        # retrive the best parameters
        max_s, max_l, max_r = np.unravel_index(
            corr_all.argmax(), corr_all.shape
        )
        
        
        self.lambda_best = self.lambdas[max_l]
        self.rsL2_best = self.rsL2s[max_r]
        self.spectrum_best = self.spectrums[max_s]

        print(f"Best parameters are: spectrum={self.spectrum_best}, lambda={self.lambda_best}, rsL2={self.rsL2_best} " )
        
        max_corr = corr_all[max_s, max_l, max_r]
        print(f"The cross validation performance is {round(max_corr, 3)} ")

        

    
        return(max_corr)
    
    def fit(self, D, P, Y):
        
        self.D = D.values.astype('double')
        self.P = P.values.astype('double')
        self.Y = Y.values.T.astype('double')
    
        
        self.gene_name = list(D.index)
        self.TF_name = list(D.columns)
        self.protein_name = list(P.columns)
        self.cell_name = list(P.index)
        
        self._cv(D, P, Y)
        
        self._fit(D, P, Y)
        

    def _fit(self, D, P, Y,lamda=None, rsL2=None,  spectrum=None ):

        """ trains a SPaRTAN model

        Parameters
        ----------
        D : array of shape (N, Q)
            The data matrix with N genes and Q TFs

        P : array of shape (M, S)
            The data matrix with M cells and S proteins

        Y : array of shape (N, M)
            The data matrix with N genes and M cells

        lamda : float > 0, default=0.001
            LASSO regularization for linear regression

        rsL2 : float > 0, default=0.001
            ridge regularization for linear regression

        corrtype: string, default='spearman'
            correlation type used to evaluate the performance

        """
        
        if lamda is None:
            lamda = self.lambda_best
        
        if rsL2 is None:
            rsL2 = self.rsL2_best
        
        if spectrum is None:
            spectrum = self.spectrum_best
        
                
        spectrumA = 1
        spectrumB = spectrum
       
        # transformation
        A = self.Y.T @ self.D
        B = self.P.T
        Y = self.Y.T @ self.Y

        # SVD(A) SVD(B)
        UA, SA, VhA = np.linalg.svd(A)
        VA = VhA.T
        UB, SB, VhB = np.linalg.svd(B)
        VB = VhB.T

        a_cum_spectrum = np.cumsum(SA) / sum(SA)
        b_cum_spectrum = np.cumsum(SB) / sum(SB)

        da = np.nonzero(a_cum_spectrum >= spectrumA)[0][0] + 1
        db = np.nonzero(b_cum_spectrum >= spectrumB)[0][0] + 1

        Ua = UA[:, :da]
        Sa = SA[:da]
        Va = VA[:, :da]

        Ub = UB[:, :db]
        Sb = SB[:db]
        Vb = VB[:, :db]

        Yv = (Y.T).flatten()

        Vb = Vb.copy(order='C')
        Ua = Ua.copy(order='C')

        L = krnP.kron(Vb, Ua)

        d = np.eye(Y.shape[0], Y.shape[1])
        cidex = np.where(d.flatten() != 0)
        diag = np.array(cidex, dtype=np.int32).flatten()

        # make it c-like contiguous array
        Yv = Yv.copy(order='C')
        diag = diag.copy(order='C')

        L, Yv = krnP.removeDiagC(L, Yv, diag)

        opts = dict()
        opts['rsL2'] = rsL2

        # reshape Yv to 2darry
        Yv = Yv.reshape(Yv.shape[0], 1)
        beta, b = leastR.LeastR(L, Yv, lamda, opts)

        del L, Yv
        gc.collect()

        self.beta = beta
        self.Ua = Ua
        self.Ub = Ub
        self.Sa = np.diag(Sa)
        self.Sb = np.diag(Sb)
        self.Va = Va
        self.Vb = Vb

    def _ar_model2w(self):
        # converts a trained model to W
        m1 = self.Va
        m2 = np.linalg.pinv(self.Sa)
        m3 = self.beta.reshape(self.Va.shape[1], self.Ub.shape[1], order="F")
        m4 = np.linalg.pinv(self.Sb)
        m5 = self.Ub.T
        ww = m1 @ m2 @ m3 @ m4 @ m5
        return ww

    def _ar_reconstruction(self, pred_test=None):
        """ reconstruction function
        Parameters
        ----------
        pred_test: prediction on test data

        """
        A = self.Y.T @ pred_test
        B = scipy.linalg.orth(self.Y)
        cm = scipy.linalg.lstsq(B, self.Y)[0]
        ct = scipy.linalg.lstsq(cm.T, A)[0]
        pred = B @ ct
        return pred

    def predict(self, P_test):
        """ predict target gene expression

        Parameters
        ----------
        P_test: Protein expression on test data

        Returns
        -------
        Y_pred: array of shape (C', N)
                The predicted Y matrix on test data set which has C' cells and N genes

        """
        
        Pv_test = P_test.values
        
        w = self._ar_model2w()
        pred = self.D @ (w @ Pv_test.T)

        aff_rec = self._ar_reconstruction(pred)

        Y_pred = aff_rec.T
        
        Y_pred = pd.DataFrame(Y_pred, index=P_test.index, columns=self.gene_name)

        return( Y_pred)
    
   

    def get_pred_score(self, P_test, Y_test, plot=False):
        """ get the correlation between predicted Y_pred and Y_test

        Parameters
        ----------
        
        Y_test: array of shape (C',N)
               gene expression test data with C' cells and N genes
        plot: whether to plot the correlation between Y_pred and Y_test, default is False


        Returns
        -------
        corr: float 0 <= value <= 1
              spearman/pearson corrlatioin between flattened self.Y_pred and Y_test

        """
        
        Y_pred = self.predict(P_test)
        
        if self.corrtype == 'spearman':
            corr = stats.spearmanr(Y_test.values.ravel(order='F'), Y_pred.values.ravel(order='F'))[0]
        else:
            corr = stats.pearsonr(Y_test.values.ravel(order='F'), Y_pred.values.ravel(order='F'))[0]

        if plot:
            plt.plot(Y_test.values.ravel(order='F'), Y_pred.values.ravel(order='F'),
                     linestyle='none', marker='+')
            plt.title('reconstruction of Y test, corr={:.2f}'.format(corr))

        return corr

    def get_W(self):
        # get coefficient matrix

        self.W = self._ar_model2w()
        return self.W

    def get_projP(self, Y=None):
        """ get projected protein expression

        Parameters
        ----------
        Y:  array of shape (optional, default is (N, M) )
            input gene expression with N genes and M cells

        Returns
        -------
        projP: array of shape (M, S)
               projected protein expression with M cells and S proteins

        """

        if Y is None:
            Yv = self.Y
        else:
            Yv = Y.values.T
            
        W = self._ar_model2w()
        
        projP = (Yv.T @ self.D @ W).T
        return (pd.DataFrame(projP, index=Y.index, columns=self.protein_name))

    def get_projTF(self, P=None):
        
        """ get projected TF activity
        Parameters
        ----------
        P: array of shape (optional, default is (M, S) )
           input protein expression with M cells and S proteins

        Returns
        -------
        projTF:  array of shape (Q, M)
            projected TF activities with Q TFs and M cells

        """
        
        
        if P is None:
            Pv = self.P
        else:
            Pv = P.values
            
        W = self._ar_model2w()
        
        projTF = (W @ Pv.T).T
        
        projTF = pd.DataFrame(projTF, index=P.index, columns=self.TF_name)
        
        return (projTF)
    
    
    def get_tf_protein_cor(self, P=None):
        '''
        Finds correlation between TFs and surface proteins
        Parameters
        ----------
        P : array_like, optional
            Protein expression matrix to use for correlations, if different from training data.

        Returns
        -------
        X : pd.DataFrame
            Correlation matrix between TFs and surface proteins.
        '''
        
        
        if P is None:
            Pv = self.P
            P = pd.DataFrame(Pv, index=self.cell_name, columns=self.protein_name)
 
        tf=self.get_projTF(P)

        # X=scipy.stats.pearsonr(
        #     tf_activity,
        #     P,
        #     axis=0)[0]



        X=scipy.stats.zscore(tf.values, axis=0).T.dot(
            scipy.stats.zscore(P, axis=0)
        )/P.shape[0]

        if self.TF_name is not None and self.protein_name is not None:
            X=pd.DataFrame(
                X,
                index=self.TF_name,
                columns=self.protein_name
            )

        return X
