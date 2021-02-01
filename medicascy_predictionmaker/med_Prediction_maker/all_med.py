import pandas as pd
# MODULE FILE
# returns df of exon list with specified MEDICASCY cutoff
# cutoff parameters (disease, treshold)
# run time: ~1 sec

def get_medicascy(disease_num, threshold):
    df = pd.read_csv("medicascy_predictionmaker/diseasefiles/{}_ef.txt".format(disease_num), sep = ' ', names = ['protein', 'num1','score1','num2', 'score2', 'score3'])
    df = df[df['score3'] >= threshold]

    return df



# print(get_medicascy(2022,30))
