import string
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

##########------------------   IMPORTANT   ------------------##########
# This script uses matplotlib and numbpy in order to generate the graphs of
# Strong scaling speedup and efficiency from the general output format that 
# was given to us in the sequentail execution of multistart_hooke code. 
# Code is not optimized but we opted to generalize it as much as possible. 
# pip3 install -r requirements.txt may be needed.
# The extracted value is the elapsed time of each file and the number of threads 
# the execution used.
# IMPORTANT: the naming of the files format must be TYPE.NUMBEROFTHREADS-CORES.txt 
# The python executable must be at the same directory as the results of our code
# the output of the sequential code must be also included as sequentialt.txt or change it below
sequential = "sequential.txt"

##########------------------     Append all txts from current Directory and extract sequential speed  ------------------##########
#Get the resulting txts from current directory
txts = glob.glob('*.txt')

#extract sequential speed and problem dimensions
seqTime = float()
ntrials = sys.argv[1]

if sequential in txts:
    with open(sequential,"r") as file:
        listoflines = file.readlines()
        seqTime = float(listoflines[1].split(" ")[len(listoflines[1].split(" ")) -2])
        file.close()
    txts.remove(sequential)
else : 
    print("No input of the sequential code found")
    exit(1)

##########------------------  Initialize ThreadCount   ------------------##########
#threadCount must be a set of the number of threads used 
threadCount = set()
for each in txts:
    threadCount.add(int(each.split(".")[1]))
threadCount = list(sorted(threadCount))

##########------------------  Initialize InputResults   ------------------##########
# inputResults will be the created dataframe
# The second line of each txt should end with #(elapsedTime) s

inputResults = []

for each in txts:
    founded = False
    partitionedEach = each.split(".")
    #out is each list of inputResults
    for out in inputResults:
        if partitionedEach[0] == out[0]:
            founded = True
            index = threadCount.index(int(partitionedEach[1]))
                        
            with open(each,"r") as file:
                listoflines = file.readlines()
                out[index+1] = seqTime/float(listoflines[1].split(" ")[len(listoflines[1].split(" ")) -2])
            file.close()
            break
    if founded is False :
        #file will be a list created for each file in txts and appended to inputResults
        inp = [0]*(len(threadCount)+1)
        inp[0]= partitionedEach[0]
        index = threadCount.index(int(partitionedEach[1]))
 
        with open(each,"r") as file:
                listoflines = file.readlines()
                inp[index+1] = seqTime/float(listoflines[1].split(" ")[len(listoflines[1].split(" ")) -2])
                file.close()
        inputResults.append(inp)

##########------------------  Print Results  ------------------##########
print(threadCount)
for inp in inputResults:
    print(inp)


def plotResults(idealDiag, inputlines, label):
    ##########------------------     Linestyles, colors and markers   ------------------##########
        #3 available line formats 
        linestyle_tuple = [
                ('loosely dashed',        (0, (5, 10))),   #for the ideal lines
                ('straight line',              'solid'),                   
                ('densely dashdotted',    (0, (3, 1, 1, 1))),          
                ('densely dashed',        (0, (5, 1))),                   
                ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1))),
            ]
        

        #8 available line colors
        mycolors = ["lightcoral","olive","seagreen","turquoise","steelblue","orchid","indigo","teal"]
        #4 available markers
        mymarkers = [ "o", "d", ">", "<"]

    ##########------------------     Create the plot    ------------------##########
        
        plt.rcParams['figure.figsize'] = [10, 14]
        plt.style.use("fivethirtyeight")
        fig, (ax1, ax2 ) = plt.subplots(2) #create a figure containing a single axes

    # Add labels-title of ax1
        ax1.set_title("Speedup for strong scaling {0}".format(label))
        ax1.set_ylabel("Speedup")
        ax1.set_xlabel("# threads-cores")
        ax1.set_ylim( 0, idealDiag[len(idealDiag)-1])
        #ax1.set_xlim( 0, idealDiag[len(idealDiag)-1])
        ax1.set_xticks(idealDiag)
        yaxis = np.linspace(idealDiag[0],idealDiag[len(idealDiag)-1],2*len(idealDiag)).round(0)
        #yaxis = idealDiag
        ax1.set_yticks(yaxis)

        #create minor ticks for axis 1
        minor_yticks = []
        y= yaxis[0]
        for x in yaxis[1:]:
            minor_yticks.append((x+y)/2)
            y = x
        ax1.set_yticks( minor_yticks, minor= True )

    # Add lables-title of ax2
        yaxis = np.linspace(0.5,1.05,2*len(idealDiag)).round(2)
        ax2.set_title("Efficiency for strong scaling {0}".format(label))
        ax2.set_ylabel("Efficiency")
        ax2.set_xlabel("# threads-cores")
        ax2.set_ylim( yaxis[0], yaxis[len(yaxis)-1])
        #ax2.set_xlim( 0 , idealDiag[len(idealDiag)-1])
        ax2.set_xticks(idealDiag)
        ax2.set_yticks(yaxis)

    # Plotting ideal lines
        ax1.plot(idealDiag, idealDiag, linewidth = 0.8, linestyle = linestyle_tuple[0][1], color= "grey", alpha=0.5 )
        ax2.plot([ idealDiag[0],idealDiag[len(idealDiag)-1]],[1,1],linewidth = 0.8, linestyle= linestyle_tuple[0][1], color= "grey", alpha=0.5 )

    ##########------------------     Plot the input     ------------------##########

        # x is a tuple where x[0] is the index and x[1] is the list
        for x in enumerate(inputlines):
            ax1.plot(idealDiag, x[1][1:], linewidth= 0.8, marker='.', markerfacecolor = mycolors[ x[0] % len(mycolors) ],label= x[1][0], linestyle = "solid", color =  mycolors[ x[0] % len(mycolors) ])
            
            efficList = []
            #y is a tuple where y[0] is the index of the list and y[1] each element of the list
            for y in enumerate(x[1][1:]):
                efficList.append( y[1] / idealDiag[y[0]])

            ax2.plot(idealDiag, efficList, linewidth= 1.5, marker= mymarkers[x[0] % len(mymarkers)], markerfacecolor = mycolors[ x[0] % len(mycolors) ], label= x[1][0], linestyle = linestyle_tuple[ (x[0]+1) % len(linestyle_tuple) ][1], color =  mycolors[ x[0] % len(mycolors) ])

        ax1.legend(loc="upper left")
        ax2.legend(loc="upper right")

        #save the image to your working directory
        plt.savefig("Speedup and Efficiency -{0} Trials -{1} Threads -{2} Aggresion Scale.png".format(sys.argv[1],sys.argv[2],sys.argv[3]))
        plt.clf()

##########------------------     Call plotResults for the input data    ------------------##########

label =  " - " + str(ntrials) + " trials"
plotResults(threadCount,inputResults, "")