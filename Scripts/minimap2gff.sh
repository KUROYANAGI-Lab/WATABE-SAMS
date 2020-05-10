#!/bin/bash
# generate gff annotations from minimap2 bam output by using Pinfish

# input
while getopts 'b:o:c:p:' optname
  do
    case $optname in
      b)
        bamFile=$OPTARG;;
      o)
        gffFile=$OPTARG;;
      c)
        clusterFile=$OPTARG;;
      p)
        pinfishDir=$OPTARG;;
    esac
  done
  

# minimap2 bam to gff
"$pinfishDir"/spliced_bam2gff/spliced_bam2gff -M -s "$bamFile" |\

# cluster gff
"$pinfishDir"/cluster_gff/cluster_gff -c 3 -d 2 -e 30 -a "$clusterFile" |\

# collapse gff
"$pinfishDir"/collapse_partials/collapse_partials -d 2 -e 30 -f 1000 > "$gffFile"
