Grit Documentation

Description
Predict transcription factor binding sites for orthologue genes using mixed Student's t-test statistics.

Homepage
http://www.thua45.cn/grit

Release
Source code     2.4.3         08/14/2021
Source code     2.4.4         12/14/2021
Binary       Windows, Mac OS, Linux

Install
Go the source folder, and run "make" command, a binary will produced in the "bin" folder

Requirement
Minimum 16GB RAM, if you install from the source code, the g++ complier is also required.

Usage
grit -m motif -i homoseq -b bgseq [-n sampling_times] [-z sampling_size] [-s rs] [-t pvalue] [- pscore] -o output

Options
-m PWMs for transcription factors
-i putative promoter sequence for orthologues genes
-b background sequences
-n sampling times, how many times sampling from background sequences, default = 10
-z sampling size, default = 200
-s raw score, 1 for RS1, 2 for RS2, default = 1
-t p-value threshold, TFBS with p-value less than will p-value threshold be reported, default = 0.05
-p p-score threshold, TFBS with p-score less than will be reported p-score threshold, default = 0
-o output, output file name

Data
motif file: Jaspar-2019.txt, Jaspar-2020.txt
promoter seq file: homoseq-v100-part1.txt, homoseq-v100-part2.txt, homoseq-v100-part3.txt, homoseq-v100-part4.txt, homoseq-v100-part5.txt, homoseq-v100-part6.txt, homoseq-v100-part7.txt, homoseq-v100-part8.txt
bakground seq file: bgseq-rdm2000.txt
results (Jaspar-2019): result-v100-part1.txt, result-v100-part2.txt, result-v100-part3.txt, result-v100-part4.txt, result-v100-part5.txt, result-v100-part6.txt, result-v100-part7.txt, result-v100-part8.txt
results (Jaspar-2020): result-v100-j20-part1.txt, result-v100-j20-part2.txt, result-v100-j20-part3.txt, result-v100-j20-part4.txt, result-v100-j20-part5.txt, result-v100-j20-part6.txt, result-v100-j20-part7.txt, result-v100-j20-part8.txt

Example
An example run should like: grit -m Jaspar-2019.txt -i homoseq-v100-part1.txt -b bgseq-rdm20000.txt -n 10 -z 200 -s 1 -t 0.05 -p 0 -o result-v100-part1.txt
This command took three input files: Jaspar-2019.txt, homoseq-v100-part1.txt, bgseq-rdm20000.txt. After finished run it will produce an output file named: result-v100-part1.txt

Assessment
Benchmark using ReMap datasets: Benchmark-ReMap.xlsx
Benchmark using public available datasets: Benchmark-ChIP-Atlas.xlsx

References

Tinghua Huang, Hong Xiao, Qi Tian, Zhen He, Min Yao. Identification of upstream transcription factor binding sites in orthologous genes using mixed Student's t-test statistics.

Contact
Dr. Tinghua Huang, thua45@126.com
Dr. Min Yao, minyao@yangtzeu.edu.cn