import pandas as pd
import iterate_diseases as ids

# grid program

tuple_list = []
for i in range(4,10):
    i = i/10
    for y in range(1,100,5):
        tuple_list.append((i,y))
print(tuple_list)
#make grid of correlations at each cutoff
df = pd.DataFrame()
for (entp,med) in tuple_list:
    correlation = ids.iterate_diseases(entp,med)
    df.loc[entp,med] = correlation
print(df)
df.to_csv("medicascy_predictionmaker/grids/grid_alltrainingTEST.txt", sep="\t")
