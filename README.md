## About DAH:
DAH is a genome assembly tool that is based on the de Bruijn graph data structure and written in the Go programming language. It is designed to use short reads to generate long contigs, and it includes features such as suggestions for selecting the optimal k-mer size and example input data. DAH is intended to be more than just a traditional assembler, as it provides additional features and resources to help users with their genome assembly projects.

## Supported data type
DAH is currently designed to work with Illumina short reads and genomes. When given short reads as input, it will first determine the appropriate k-mer length and then use that information to construct a simplified graph. From this graph, it will output contigs in a fasta file. When given a genome as input, DAH will generate pseudo reads based on user-specified coverage levels and use these synthetic reads to perform assembly. Overall, the goal of DAH is to generate high-quality contigs from either short reads or a genome using a de Bruijn graph-based approach.

## DAH pipeline
DAH consists of several components or modules, including:

` K-mer length selection: This module allows users to choose the k-mer length to use in the assembly process. There are two modes of input: default, in which DAH suggests the optimal k-mer length, or range, in which the user specifies a range of k-mer lengths to consider.

` Graph simplification and error correction: This module is responsible for processing the input data and constructing a de Bruijn graph, as well as simplifying the graph and correcting any errors that may be present.

` Visualization: This module provides users with a plot of the k-mer frequency distribution, which can be useful for understanding the properties of the input data and the resulting assembly.

The combination of these modules enables users to generate contigs from short reads or a genome using a de Bruijn graph-based method.
## External resorces needed for DAH
  To visualize the results of DAH, we used a tool called [Bandage](https://rrwick.github.io/Bandage/). Bandage is a tool for visualizing graph files using the GFA format, which is a file format that is used to represent graphs in genome assembly. Three GFA files will be produced following assembly, one for each step of the process: building the de Bruijn graph, performing the chains merging, and clipping tips. These GFA files can be visualized using Bandage to understand and analyze the assembly process and its results.
 
## DAH sample input and command line examples:
DAH includes a set of sample input files in the SampleInput folder. These files can be used for testing or as examples of the types of input that DAH can process. To use these sample inputs, you can use the appropriate command line options to specify the input file and any other necessary parameters. Here are some examples:
 - If the input is a reference Genome:
``` 
go run ./ genome SampleInput/genome 10 default
go run ./ genome SampleInput/genome 10 range 21 30
```
- If the input is reads:
``` 
go run ./ reads GeneratedReads.fasta 10 default
go run ./ reads GeneratedReads.fasta 10 range 21 30
```
There is also a link to the representation video:
[Link](https://drive.google.com/file/d/1bjF9t_3Yt4JLx0EOyWDgOFlTzk7ilc1j/view?usp=sharing)
