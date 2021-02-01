import pandas as pd
import make_prediction as mp
# takes about 90 sec
# from sklearn.metrics import matthews_corrcoef, confusion_matrix

oprev_list = [0.07108, 0.21481, 0.25424, 0.10592]
list = [2227,91,2646,1984]  #ovarian(2227), lymphoma(91), thyroid(2646), CF(1984), Asthma* (2350), bronchitis*(2978), colon cancer*(2168), breast cancer (2022),kidney cancer(2288), nephritis (1042)

def iterate_diseases(entp,med):
    pprev_list = []
    for dis in list:
        (disease,pprev,path,med) = mp.make_prediction(dis, entp, med) #dis, pred, pathogenicity, medicascy threshold
        pprev_list.append(pprev)
    df = pd.DataFrame()
    df['observed'] = oprev_list
    df['predicted'] = pprev_list

    df = df.corr(method='pearson')

    pprev_list = []

    return df.loc['observed', 'predicted']  #return pearson correlation


# print(iterate_diseases(0.9,0))
# print("this separates ")
# print(iterate_diseases(0.5,0.1))

# To do:
# Organize list of diseases and file numbers to run MEDICASCY on
# Get the sweet spot cutoff, and then test on the testing set.
# Get the correlation between cutoff predictions and observed and report that.
# Do the same for Know-GENE.
