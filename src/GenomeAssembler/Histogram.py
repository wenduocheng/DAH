
#This function Histogram is used to generate the kmer distribution plot by reading a text file that each element is SEPERATE BY !!!!!!![,](might change latter) as input, and draw the histoplot based on that data. filedir is the path and the name of the file that I want to draw based on

def Histogram(filedir):
    from matplotlib import pyplot as plt
    import numpy as np  
#first, read the data by the given filedir
    text_file = open(filedir,"r")
    data = text_file.read()
    data_split = data.replace('\n', ' ').split( ',' )
    data_split = data_split[:-1]
    #make sure I got the line
    print(data_split)
    text_file.close()
    #convert all the string elements into int
    x = [int(i) for i in data_split]
    #verify and print again
    print(x)
    #data set:
    #x = np.array([3,3,9,9,18,18,18,18,18,18,15,12,6,6,12,15])
    #set bin size
    bin_size = np.unique(x)

   
    plt.xlabel("frequency")
    plt.ylabel("number of distinct kmers")
    plt.xlim
    plt.hist(x, bins = bin_size, edgecolor = "black")
    plt.savefig("25.png")
    plt.show()
    
Histogram("25.txt")
  
def Histogram_xy(filedir):
    from matplotlib import pyplot as plt
    import numpy as np
    x_text_file = open(filedir,"r")
    x_data = x_text_file.read()
    x_data_split =x_data.replace('\n', ' ').split( ',' )
    x_data_split = x_data_split[:-1]
    x = [ idx for idx in range(len( x_data_split))]
    y=[int(i) for i in x_data_split]
    # bin_size = np.unique(x)
   
    plt.xlabel("Position in chromosome")
    plt.ylabel("coverage")
    plt.bar(x,y, edgecolor = "blue")
    plt.savefig("read for plot")
    plt.show()

Histogram_xy("reads for plot.txt")




    

