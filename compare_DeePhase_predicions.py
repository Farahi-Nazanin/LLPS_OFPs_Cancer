#
#
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

def read_results(infilename, column):
    mylist = []
    with open(infilename) as infile:
        for line in infile:
            if line.strip().split('\t')[column] not in ['', 'LLPS PROBABILITY']:
                if line.strip().split('\t')[column][1:2] != '.':
                    print(line)
                else:
                    mylist.append( float(line.strip().split('\t')[column]) )
    return mylist


scaffolds = read_results('../LLPS_prediction_drivers/DeePhase_predictions/results.txt', 1)
fusions = read_results('../LLPS_prediction_fusions/DeePhase_predictions/results_COSMIC_fusions.txt', 1)
COSMIC_census = read_results('DeePhase_predictions/results_COSMIC_census.txt', 1)
proteome = read_results('DeePhase_predictions/results.txt', 1)

names = ["Scaffolds", "Fusions", "COSMIC", "Proteome"]
lists = [scaffolds, fusions, COSMIC_census, proteome]
for i in range(0,3):
    for j in range(i+1, 4):
        U, p = mannwhitneyu(lists[i], lists[j] )
        print("Avg("+names[i]+") = "+str(np.mean(lists[i]))+ \
              "; Avg("+names[j]+") = "+str(np.mean(lists[j])) )
        print("p("+names[i]+", "+names[j]+") =", p)
        print()

meanlineprops = dict(linestyle='dotted', linewidth=1.6, color='black')
bplot = plt.boxplot(x = [scaffolds, fusions, COSMIC_census, proteome], \
            labels = ['LLPS scaffolds', 'Oncogenic fusions', 'COSMIC Census', 'Proteome'], \
            notch=True, showmeans=True, meanline=True, meanprops=meanlineprops, patch_artist=True)
# fill with colors
colors = ['firebrick', 'lightpink', 'royalblue', 'lightsteelblue']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color('black')
plt.title('Distributions of LLPS propensities')
plt.ylabel('DeePhase prediction score')
plt.axhline(0.5, c='r', linestyle='--')
plt.ylim(0,1)
plt.savefig("DeePhase_predictions/DeePhase_predictions_final.png", dpi=600)

