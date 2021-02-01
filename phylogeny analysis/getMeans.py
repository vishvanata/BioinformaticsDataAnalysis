from getScoreList import getScoreList
from idiosync import idiosync
import matplotlib.pyplot as plt


# testScores = [0.45,0.56,0.57,0.58,0.65,0.78]
meanlist = []
x = [1,2,3,4,5,6,7,8,9]
x.reverse()

for i in range(1,10):
    input = i;
    scoreList = getScoreList(idiosync(input))
    meanlist.append(sum(scoreList,0)/len(scoreList))
meanlist.reverse()
plt.scatter(x,meanlist, c ="blue");
plt.title('Scatterplot of means'.format(input))
plt.xlabel('Common to n Genes')
plt.ylabel('Mean Pathogenicity')
plt.show()


    # bins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]

    # plt.hist(scoreList, bins, histtype='stepfilled', rwidth=0.8, color='green')
    # plt.xlabel('Pathogenicity')
    # plt.ylabel('Frequency')
    # plt.title('Genes common to {} cancers'.format(input))
    # plt.savefig('unique_distributions/uniqdistr_{}.png'.format(input))
    # plt.show()
