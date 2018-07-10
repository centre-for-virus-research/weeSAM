#!/usr/local/bin/python3 -u
#'''
#wee_pysam is a python (pysam) re-write of weeSAM written by Joseph Hughes.
#
#His original script can be found here:
#https://github.com/josephhughes/Sequence-manipulation/blob/master/weeSAMv1.4
#
#This script produces the following information from an input bam or sam file:
#    -   Number of mapped reads
#    -   Average coverage
#    -   Min/Max coverage
#    -   Breadth/Depth of coverage
#along with more information.
#'''

# Import packages
import argparse
import pysam
import os
import csv
import statistics
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt



#'''
# Make argparse variables for:
#    -    sam file
#    -    bam file
#    -    cut off
#    -    out file
#    -    plot
#'''
parser = argparse.ArgumentParser(description="A pysam script to produce information from an "+
            "input BAM or SAM file.")
parser.add_argument("--sam", help="The input SAM file.")
parser.add_argument("--bam", help="The input BAM file.")
parser.add_argument("--cutoff", help="Cut off value for number of mapped reads. [Default = 0]")
parser.add_argument("--out", help="Output file name.", required=True)
parser.add_argument("--plot", help="An output file name for the coverage plot. (.pdf)")
args = parser.parse_args()

#'''
# Check if a sam or bam has been given to the command line. Only one can be given.
#'''
# Check if none exist
if args.sam is None and args.bam is None:
   print("You haven't specified a bam or a sam file. Please do this with either the"
         "\n\t--sam `file.sam`\n\t--bam `file.bam`")
   exit(1)
# Check if both exists
elif args.sam and args.bam:
    print("You've specified both a sam and bam file, you can only specify one."
        " Please use either\n\t--sam `file.sam`\n\tOR\n\t--bam`file.bam`")
    exit(1)

#'''
#Check if the user has specified a cutoff value. if not set it to zero.
#'''
if args.cutoff is None:
    args.cutoff = 0
else:
    args.cutoff = int(args.cutoff)

#'''
#A function which plots a line graph based on an input list and a title. The funciton returns the plot.
#'''
def my_plot(lst, name):
    fig = plt.figure(figsize=(8,6))
    plt.plot(lst)
    plt.xlabel("Position")
    plt.ylabel("Depth of coverage")
    plt.title(name)

    return fig

#'''
# If a user has given a sam file sort it and call it their input minus .sam plus .bam
# after this args.bam is always going to exists, so I can write code for args.bam and
# it will always work.
#'''
if args.sam:
    print("Sorting SAM file...")
    pysam.sort("-o", args.sam.rsplit(".")[0]+".bam", args.sam)
    args.bam = args.sam.rsplit(".")[0]+".bam"
#    os.remove(args.sam)

#'''
# Work on the bam file generated from sam2bam or the one given on the command line.
#'''
if args.bam:
    # Index bam file
    pysam.index(args.bam)
    #'''
    #samtools idxstats the bam file and write the results to tmp_mapped.txt.
    #Then pull out mapped unmapped reads from this file.
    #'''
    print(pysam.idxstats(args.bam), end="", file=open("tmp_mapped.txt","w"))
    with open("tmp_mapped.txt", "r") as f:
        # Dictionary which stores the values to print to output file later.
        # { NAME : REFLEN\tMAPPEDREADS }
        stat_dict = {}
        contents = csv.reader(f, delimiter="\t")
        # Counter for tot unmapped / mapped
        total_unmapped = 0
        total_mapped = 0
        for line in contents:
            total_unmapped += int(line[3])
            total_mapped += int(line[2])
            if int(line[2]) > int(args.cutoff):
                key = stat_dict.setdefault(line[0],[])
                key.append(line[1])
                key.append(line[2])

        f.close()

    #Remove the temp samtools idxstats file
    os.remove("tmp_mapped.txt")
    # Print summary stats of reads to screen.
    print("Total # Mapped Reads:\t"+str(total_mapped)+"\nTotal # Unmapped Reads:\t"+str(total_unmapped))
    print("Total # Reads:\t"+str(int(total_mapped)+int(total_unmapped)))

    #'''
    # samtools depth the bam file for coverage statistics.
    #'''
    print(pysam.depth("-d", "1000000", args.bam), end="", file=open("tmp_depth.txt", "w"))

    #'''
    # Open the file read the lines, store as a list close the file.
    #                   0   5   /   0   7
    # ! !   Ideally i would want to do this all whilst opening file and storing nothing in mem ! !
    #                   0   6   /   0   7
    # !!    Script now reads it into one dict instead of five or six
    #'''
    min_dict = {}
    with open("tmp_depth.txt", "r") as f:
        data = csv.reader(f, delimiter="\t")
        for line in data:
            key = min_dict.setdefault(line[0],[])
            key.append(int(line[2]))
        f.close()


    # Remove the temp file
    os.remove("tmp_depth.txt")

    #'''
    # Store the average depth as a dictionary. This will be needed to calculate site with 0.2 of mean cov
    # etc. Looping through the temp depth data stored in mem to see these sites.
    #'''
    avg_dict = {}
    for i in stat_dict:
        avg_dict[i] = statistics.mean(min_dict[i])

    #'''
    #Use the avg dict from above to loop through each element in dict and check whether its value is above the
    #specified threshold.
    #if it is append the value 1 to a dict list, one for each threshold.
    #Maths can then be applied to this new list to work out percentage of site above the threshold.
    #'''
    two_dict = {}
    five_dict = {}
    eighteen_dict = {}
    one_dict = {}
    percent_dict = {}
    for i in stat_dict:
        for j in min_dict[i]:
            # Get the total number of sites into a dict. You can then compare how many elements in this
            # dict against the ref len to get percent coverage
            percent_key = percent_dict.setdefault(i, [])
            percent_key.append(1)
            # If the criteria is met make the dict value a list, and append 1 to the dict value
            if j > float(avg_dict[i] * 1.8):
                key = eighteen_dict.setdefault(i, [])
                key.append(1)
            # If it isnt met create the dict value but append 0 to it. This is needed as if the criteria
            # isnt met the dict key/value isnt created thus it can't be printer later.
            else:
                key = eighteen_dict.setdefault(i, [])
                key.append(0)
            if j > float(avg_dict[i]):
                key_1 = one_dict.setdefault(i, [])
                key_1.append(1)
            else:
                key_1 = one_dict.setdefault(i, [])
                key_1.append(0)
            if j > float(avg_dict[i] * 0.5):
                key_5 = five_dict.setdefault(i, [])
                key_5.append(1)
            else:
                key_5 = five_dict.setdefault(i, [])
                key_5.append(0)
            if j > float(avg_dict[i] * 0.2):
                key_2 = two_dict.setdefault(i, [])
                key_2.append(1)
            else:
                key_2 = two_dict.setdefault(i, [])
                key_2.append(0)

    #'''
    #Open a file, called whatever the user specified and print to it.
    #Open a html file if args.plot is given to argparse
    #'''
    out_file = open(args.out, "w")
    if args.out:
        html_file = open(args.out.split(".")[0]+".html","w")
        html_list = []

    # Print the file headers
    print("Ref_Name\tRef_Len\tMapped_Reads\tBreadth\t%_Covered\tMin_Depth\tMax_Depth\tAvg_Depth\t"
          "Std_Dev\tAbove_0.2_Depth\tAbove_0.5_Depth\tAbove_1_Depth\tAbove_1.8_Depth\tVariation_Coefficient", file=out_file)
    for i in stat_dict:
        print("Processing:\t"+str(i)+"\t....")
        if args.plot:
            fig = my_plot(min_dict[i], i)
            fig.savefig(str(i) + ".svg", format="svg")
            html_list.append("""
            <tr>
                <td><a href={image}>{ref_name}</a></td>
                <td>{ref_len}</td>
                <td>{mapped_reads}</td>
                <td>{breadth}</td>
                <td>{covered}</td>
                <td>{min}</td>
                <td>{max}</td>
                <td>{avg}</td>
                <td>{stdev}</td>
                <td>{point2}</td>
                <td>{point5}</td>
                <td>{one}</td>
                <td>{one_eight}</td>
                <td>{var_co}</td>
            </tr>
            """.format(image=str(i)+".svg", ref_name=str(i), ref_len=str(stat_dict[i][0]),
                       mapped_reads=str(stat_dict[i][1]),
                       breadth=str(sum(percent_dict[i])),
                       covered=str(sum(percent_dict[i])/int(stat_dict[i][0])*100),
                       min=str(min(min_dict[i])), max=str(max(min_dict[i])),
                       avg="%.2f"%(statistics.mean(min_dict[i])),
                       stdev="%.2f"%(statistics.stdev(min_dict[i])),
                       point2="%.2f"%(sum(two_dict[i])/int(stat_dict[i][0])*100),
                       point5="%.2f"%(sum(five_dict[i])/int(stat_dict[i][0])*100),
                       one="%.2f"%(sum(one_dict[i])/int(stat_dict[i][0])*100),
                       one_eight="%.2f"%(sum(eighteen_dict[i])/int(stat_dict[i][0])*100),
                       var_co="%.2f"%(statistics.stdev(min_dict[i])/statistics.mean(min_dict[i]))))

        if args.out:
            # Print contig / seq name   ref len     reads mapped
            print(str(i)+"\t"+str(stat_dict[i][0])+"\t"+str(stat_dict[i][1])+"\t"+
            # Print the breadth and the percentage covered
            str(sum(percent_dict[i]))+"\t"+str(int(sum(percent_dict[i])/int(stat_dict[i][0])*100))
            # Print the minimum and maximum depth values for each seq / contig
            +"\t"+str(min(min_dict[i]))+"\t"+str(max(min_dict[i]))+
            # Print the mean value of depth for each seq / contig to two dp and the std devitaton
            "\t"+"%.2f"%(statistics.mean(min_dict[i]))+"\t"+"%.2f"%(statistics.stdev(min_dict[i]))+"\t"+
            # Print the % of sites which meet the criteria specified to two dp
            #       sites above 0.2 mean coverage
            "%.2f"%(sum(two_dict[i])/int(stat_dict[i][0])*100)+"\t"+
            #       sites above 0.5 mean coverage
            "%.2f"%(sum(five_dict[i])/int(stat_dict[i][0])*100)+"\t"+
            #       sites above mean coverage
            "%.2f"%(sum(one_dict[i])/int(stat_dict[i][0])*100)+"\t"+
            #       sites above 1.8 mean coverage.
            "%.2f"%(sum(eighteen_dict[i])/int(stat_dict[i][0])*100)+"\t"+
            # Print the variation coefficient
            "%.2f"%(statistics.stdev(min_dict[i])/statistics.mean(min_dict[i])),file=out_file)

html_str = """
<!doctype html>
<html>
<head>
<style>
table {{
    font-family: arial;
    border-collapse: collapse;
    width: 100%;
}}

td, th {{
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
}}

tr:nth-child(even) {{
    background-color: #dddddd;
}}

</style>
</head>
<body>

<h2>{title}</h2>     

<table>
    <tr>
        <th>Ref_Name</th>
        <th>Ref_Len</th>
        <th>Mapped_Reads</th>
        <th>Breadth</th>
        <th>%_Covered</th>
        <th>Min_Depth</th>
        <th>Max_Depth</th>
        <th>Avg_Depth</th>
        <th>Std_Dev</th>
        <th>Above_0.2_Depth</th>
        <th>Above_0.5_Depth</th>
        <th>Above_1_Depth</th>
        <th>Above_1.8_Depth</th>
        <th>Variation_Coefficient</th>
    </tr> 

"""

html_str = html_str.format(title="weeSAM output for file:\t"+str(args.bam))
html_file.write(html_str)
for i in html_list:
    print(i, file=html_file)
html_file.write("\t</table>\n</body>")
