### _wee\_pysam_
```  
              _)_
           .-'(/ '-.
          /    `    \
         /  -     -  \
        (`  a     a  `)
         \     ^     /
          '. '---' .'
          .-`'---'`-.
         /           \
        /  / '   ' \  \
      _/  /| S A M |\  \_
     `/|\` |+++++++|`/|\`
          /\       /\
          | `-._.-` |
          \   / \   /
          |_ |   | _|
          | _|   |_ |
          (ooO   Ooo)


``` 

### _usage:_ ### 
``wee_sam.py { --bam file.bam OR --sam file.sam } --out output.txt `` 

### _compiled zip usage:_ ###   
_the .zip contains all the required modules so should run on any system with python3 installed_  
``python3 weesam.zip { --bam file.bam OR --sam file.sam } --out output.txt ``   

_command line options:_   
`--bam` : The input bam file.  
`--sam` : The input sam file.  
`--out` : The name of your output file.  
`-h`    : Help  

_File output fields:_  
`Ref_Name` : The identifier of the reference.  
`Ref_Len` : The length in bases of each reference.  
`Mapped_Reads` : Number of reads mapped to each reference.  
`Breadth` : The number of sites in the genome covered by reads.  
`%_Covered` : The percent of sites in the genome which have coverage.  
`Min_Depth` : Minimum read depth observed.  
`Max_Depth` : Max read depth observed.  
`Avg_Depth` : Mean read depth observed.  
`Std_Dev` : Standard deviation of the mean (Avg_Depth).  
`Above_0.2_Depth` : Percentage of sites which have greater than 0.2 * Avg_Depth.  
`Above_0.5_Depth` : Percentage of sites which have greater than 0.5 * Avg_Depth.  
`Above_1_Depth` : Percentage of sites which are above Avg_Depth.  
`Above_1.8_Depth` : Percentage of sites which have greater than 1.8 * Avg_Depth.  
`Variation_Coefficient` : The mean of Std_Dev of the mean.  


