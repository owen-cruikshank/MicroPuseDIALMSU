
import decimal
import numpy as np
import pickle

def TIPS_2017_Owen(mol,iso,T):
    niso = [0,9,13,18,5,9,4,6,3,2,1,2,2,3,2,4,4,2,2,5,3,2,3,3,2,1,3,2,1,
    2,1,3,1,1,1,2,1,2,2,1,2,4,1,1,6,2,4,1,2,2,3,1,1,4]

    mol = str(mol)
    iso = str(iso)
    file = 'QTpy/'+mol+'_'+iso+'.QTpy'

    QTdict = {}
    with open(file, 'rb') as handle:
        QTdict = pickle.loads(handle.read())

    QT = T
    with np.nditer(QT, op_flags=['readwrite']) as QTT:
        for x in QTT:
            if(x==int(x)):
                key = str(int(x))
                x[...] = float(QTdict[key])
            else:
                key=str(int(x))
                Q1 = float(QTdict[key])
                key=str(int(x+1))
                Q2 = float(QTdict[key])
                x[...] = Q1+(Q2-Q1)*(x-int(x))

    return(QT)

if __name__ == '__main__':

    QT = TIPS_2017_Owen(mol,iso,T)
