"""
class of SPaRTAN model.

This script requires Pandas, numpy, Scipy, matplotlib and scikit-learn to be installed within the running environment

This script requires Cython modules present in the current directory
"""

import numpy as np
import pandas as pd
import cythKronPlus as krnP
import cythLeastR as leastR
import scipy.linalg
import gc
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.model_selection import KFold

    
class pySPaRTAN:
    """
    The major class for SPaRTAN, establishing an interaction matrix between
    surface proteins (P) and TFs (D) that predicts target gene expression (Y).

    Methods
    -------
    __init__(self, lambdas, rsL2s, spectrums, corrteyp, n_fold)
        Initiate the class object 
    
    fit(self, D, P, Y)
        Cross validates and train a SPaRTAN model

    _cv(self, D, P, Y)
         Determine the optimal hyperparameters with cross-validation and grid search
    
    _fit(self, D, P, Y, lamda=None, rsL2=None,  spectrum=None )
         train a SPaRTAN model with specific hyperparameters
     
         
    _ar_model2w(self):
        Convert trained model to interaction matrix W

    _ar_reconstruction(self, pred_test=None):
        Reconstruct the prediction

    predict(self, P_test=None):
        Predict gene expression given protein expression.

    get_pred_score(self, Y_pred, Y_test, plot=False):
        Compute the correlation between predicted Y_pred and Y_test

    get_W(self):
        Extract coefficient matrix W

    get_projP(self, Y=None):
        Obtain projected protein expression

    get_projTF(self, P=None):
        Obtain projected TF activity
        
    get_tf_protein_cor(self, P=None)
        Compute pair-wise correlations between inferred TF activity and protein expression
    """
    
    def __init__(self, lambdas=None, rsL2s=None,  spectrums=None, corrtype='pearson', n_fold=5):
        
        """
        Function:
        ----------
        pySPaRTAN initiation
        
        Parameters:
        ----------
        lambdas: list of float numbers
            hyperparameter of pySPaRTAN to tune
        rsL2s: list of float numbers
            hyperparameter of pySPaRTAN to tune
        spectrums: list of float numbers, default=[0.7]
            hyperparameter of pySPaRTAN to tune, needs less tune, fixed to 0.7
        corrtype: string
            type of correlation method used to calculate correlation between predicted and true gene expression
        n_fold: integer
            cross-validation fold
            
        Returns
        -------
        None.
        """

        self.rsL2s = rsL2s
        self.lambdas = lambdas
        self.spectrums = spectrums
  
        self.corrtype = corrtype
        self.n_fold=n_fold
       
    def fit(self, D, P, Y):
        
        """
        Function:
        ----------
        User interface for training the model
        
        Parameters:
        ----------
        D: Pandas dataframe
            gene-tf regulation with gene X TF 
        P: Pandas dataframe
            Protein expression with cell X protein
        Y: Pandas dataframe
            Gene expression with cell X gene
            
            
        Returns
        -------
        None.
        """
    
        self._cv(D, P, Y)
        
        self._fit(D, P, Y)
            
       
    def  _cv(self, D, P, Y):
        
        """
        Function:
        ----------
        Determine the optimal hyperparameters with cross-validation and grid search. Only used by model itself.
        
        Parameters:
        ----------
        D: Pandas dataframe
            gene-tf regulation with gene X TF 
        P: Pandas dataframe
            Protein expression with cell X protein
        Y: Pandas dataframe
            Gene expression with cell X gene
            
        Returns
        -------
        None.
        """

        lensps = len(self.spectrums)
        lenlambdas = len(self.lambdas)
        lenrsL2s = len(self.rsL2s)
  
          
        corr_all = np.zeros((lensps, lenlambdas, lenrsL2s))
        
        # for a in range(0, lenspAs):
        for s in range(lensps) :
            for l in range(lenlambdas):
                for r in range(lenrsL2s):
                    
                    # cross-validation
                    print(
                        "cross validating  spectrum={}, lambda={}, rsL2={}".format(
                            self.spectrums[s], self.lambdas[l] ,self.rsL2s[r]
                        )
                    )

                    Y_pred_all = Y * 0
                    
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
    
    

    def _fit(self, D, P, Y,lamda=None, rsL2=None,  spectrum=None ):
        
        """
        Function:
        ----------
        Train pySPaRTAN model with specific data and hyperparameters. Only used by model itself.
        
        Parameters:
        ----------
        D: Pandas dataframe
            gene-tf regulation with gene X TF 
        P: Pandas dataframe
            Protein expression with cell X protein
        Y: Pandas dataframe
            Gene expression with cell X gene
        lamda: float
            hyperparameter passed in to train the model
        rsL2: float
            hyperparameter passed in to train the model
        spectrum: float
            hyperparameter passed in to train the model
            
        Returns
        -------
        None.
        """

        self.D = D.values.astype('double')
        self.P = P.values.astype('double')
        self.Y = Y.values.T.astype('double')
    
        
        self.gene_name = list(D.index)
        self.TF_name = list(D.columns)
        self.protein_name = list(P.columns)
        self.cell_name = list(P.index)
        
        
        if lamda is None:
            lamda = self.lambda_best
        
        if rsL2 is None:
            rsL2 = self.rsL2_best
        
        if spectrum is None:
            spectrum = self.spectrum_best
        
                
        spectrumA = 1
        spectrumB = spectrum
        
        Dv = self.D
        Pv = self.P
        Yv = self.Y
       
        # transformation
        A = Yv.T @ Dv
        B = Pv.T
        Y = Yv.T @ Yv

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

        Yt = (Y.T).flatten()

        Vb = Vb.copy(order='C')
        Ua = Ua.copy(order='C')

        L = krnP.kron(Vb, Ua)

        d = np.eye(Y.shape[0], Y.shape[1])
        cidex = np.where(d.flatten() != 0)
        diag = np.array(cidex, dtype=np.int32).flatten()

        # make it c-like contiguous array
        Yt = Yt.copy(order='C')
        diag = diag.copy(order='C')

        
        L, Yt = krnP.removeDiagC(L, Yt, diag)

        opts = dict()
        opts['rsL2'] = rsL2

        # reshape Yt to 2darry
        Yt = Yt.reshape(Yt.shape[0], 1)
        beta, b = leastR.LeastR(L, Yt, lamda, opts)

        del L, Yt
        gc.collect()

        self.beta = beta
        self.Ua = Ua
        self.Ub = Ub
        self.Sa = np.diag(Sa)
        self.Sb = np.diag(Sb)
        self.Va = Va
        self.Vb = Vb

    def _ar_model2w(self):
        
        """
        Function:
        ----------
        Convert a trained model to W. Only used by model itself.
        
        Parameters:
        ----------
        None
            
        Returns
        -------
        ww: 2D array with the format as TF X protein
            coefficient matrix between transcription factor and protein expression
        """
        
        m1 = self.Va
        m2 = np.linalg.pinv(self.Sa)
        m3 = self.beta.reshape(self.Va.shape[1], self.Ub.shape[1], order="F")
        m4 = np.linalg.pinv(self.Sb)
        m5 = self.Ub.T
        ww = m1 @ m2 @ m3 @ m4 @ m5
        return ww

    def _ar_reconstruction(self, pred_test=None):
        
        """
        Function:
        ----------
        Reconstruct the prediction. Only used by model itself
        
        Parameters:
        ----------
        pred_test: 2D array with the format as cell X gene
            predicted gene expression before reconstruction
            
        Returns
        -------
        pred: 2D array with the format as gene X cell
            Reconstructed gene expression which the final state of predicted gene expression
        """
        
        A = self.Y.T @ pred_test
        B = scipy.linalg.orth(self.Y)
        cm = scipy.linalg.lstsq(B, self.Y)[0]
        ct = scipy.linalg.lstsq(cm.T, A)[0]
        pred = B @ ct
        return pred

    def predict(self, P_test):

        """
        Function:
        ----------
        Predict gene expression given protein expression.
        
        Parameters:
        ----------
        P_test: Pandas dataframe with the format as cell X protein
            Input protein expression to test the prediction of gene expression
            
        Returns
        -------
        Y_pred: Pandas dataframe with the format as cell X gene
            Predicted gene expression
        """
        
        Pv_test = P_test.values
        
        w = self._ar_model2w()
        pred = self.D @ (w @ Pv_test.T)

        aff_rec = self._ar_reconstruction(pred)

        Y_pred = aff_rec.T
        
        Y_pred = pd.DataFrame(Y_pred, index=P_test.index, columns=self.gene_name)

        return( Y_pred)
    
   

    def get_pred_score(self, P_test, Y_test, plot=False):

        """
        Function:
        ----------
        Compute the correlation between predicted and true gene expression with the method saved in self.corrtype
        
        Parameters:
        ----------
        P_test: Pandas dataframe with the format as cell X protein
            Input protein expression to test the prediction of gene expression
        Y_test: Pandas dataframe with the format as cell X gene
            True gene expression 
            
        Returns
        -------
        corr: float
            correlation between predicted and true gene expression
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
        
        """
        Function:
        ----------
        Get coefficient dataframe between protein expression and transcription factor
        
        Parameters:
        ----------
        None
        
        Returns
        -------
        self.W: Pandas dataframe with the format as TF X protein
             coefficient matrix between transcription factor and protein expression   
        """

        self.W = self._ar_model2w()
        W = pd.DataFrame(self.W, index=self.TF_name, columns=self.protein_name)
        
        return(W)

    def get_projP(self, Y):

        """
        Function:
        ----------
        Obtains the predicted protein expression
        
        Parameters:
        ----------
        Y: Pandas dataframe with the format as cell X gene
            Input gene expression to obtain predicted protein expression
        
        Returns
        -------
        projP: Pandas dataframe with the format as cell X protein
             predicted protein expression
        """
     
        if Y is None:
            Yv = self.Y
        else:
            Yv = Y.values.T
            
        W = self._ar_model2w()
        
        projP = (Yv.T @ self.D @ W).T
        projP = pd.DataFrame(projP, index=Y.index, columns=self.protein_name)
        
        return (projP)

    def get_projTF(self, P):
        
        """
        Function:
        ----------
        Obtains the predicted transcription factor activity
        
        Parameters:
        ----------
        P: Pandas dataframe with the format as cell X protein
            Input protein expression to obtain predicted TF
        
        Returns
        -------
        projTF: Pandas dataframe with the format as cell X TF
             predicted transcription factor activity
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

        """
        Function:
        ----------
        Calculate pair-wise correlation(similarity) between predicted TF activity and protein expression
        
        Parameters:
        ----------
        P: Pandas dataframe with the format as cell X protein
            Input protein expression to obtain predicted TF
        
        Returns
        -------
        score: Pandas dataframe with the format as TF X protein
             correlation (similarity) between predicted TF activity and protein expression
        """
        
        if P is None:
            Pv = self.P
            P = pd.DataFrame(Pv, index=self.cell_name, columns=self.protein_name)
 
        tf=self.get_projTF(P)

 
        X=scipy.stats.zscore(tf.values, axis=0).T.dot(
            scipy.stats.zscore(P, axis=0)
        )/P.shape[0]

        if self.TF_name is not None and self.protein_name is not None:
            score=pd.DataFrame(
                X,
                index=self.TF_name,
                columns=self.protein_name
            )

        return(score)
