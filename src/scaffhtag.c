/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of scaff10x pipeline.                                 *
 *                                                                          *
 *  Scaff10x is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 400 
#define Max_N_NameBase2 400 
#define Max_N_Pair 100
static char bindir[2000];
static char tmpdir[2000];
static char **S_Name;
static int *insert_siz,*insert_dev,*core_list,*ctg_list,*ctg_head,*read2contig;
static int *readIndex;

/* SSAS default parameters   */
static int n_group=0;
static int num_reads=10;
static int num1_reads=8;
static int num2_reads=10;
static int num_links=8;
static int num1_links=8;
static int num2_links=10;
static int len_block=50000;
static int len_edges=50000;
static int len_matrx=2000;
static int min_edge = 5;
static int file_tag = 0;
static int gap_len = 100;
static int uplinks = 50;
static int mscore = 20;
static int n_longread = 1;
static int run_align = 1;
static int plot_10x = 0;
static int mkdup_flag = 1;
static int min_len = 3000;
static int sam_flag = 0;
static int rc_flag = 0;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

void
RunSystemCommand(char *cmd)
{
    int ret;
    if ((ret = system(cmd)) != 0)
    {
        fprintf(stderr, "Error running command: %s\n", cmd);
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char **argv)
{
    int i,nSeq,args;
    char *st,*ed;
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    void ArraySort_String2(int n,char **Pair_Name,int *brr);
    fasta *seq;
    FILE *fp,*namef,*namef2;
    int size_file;
    int m_score,n_nodes,n_debug,num_sigma;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Read_Index(int nSeq,char *namefile);
    void Read_Group(fasta *seq,int nSeq,int nRead,int cindex);
    void File_Output(int aaa);
    void Memory_Allocate(int arr);
    char tempa[2000],tempc[2000],syscmd[2000],workdir[2000];
    char file_tarseq[2000],file_scaff[2000],file_sfagp[2000],file_datas[2000],file_plot10x[2000],file_cover[2000],file_falgn[2000];
    char file_read1[2000],file_read2[2000],samname[500],bamname[500],toolname[500],datname[500],lenplot10x[200],mkdup[200];
    int systemRet = system (syscmd);
    int systemChd = chdir(tmpdir);
    pid_t pid;

    seq=NULL;
    
    if(argc < 2)
    {
         printf("Program: scaff10x - Genome Scaffolding using 10X Chromium Data\n");
         printf("Version: 4.2\n");
         printf("\n");
         
         printf("Usage: %s -nodes 30 -matrix 2000 -reads 4 -score 20 -edge 50000 -link 3 -block 50000 -mkdup Dupmarked.bam -plot barcode_lengtg.png <input_assembly_fasta/q_file> <Input_read_1>> <Input_read_2> <Output_scaffold_file>\n",argv[0]);
         printf("       nodes    (30)    - number of CPUs requested\n");
         printf("       matrix   (2000)  - relation matrix size\n");
         printf("       reads    (10)    - step 1 and 2: minimum number of reads per barcode\n");
         printf("       link     (8)     - step 1 and 2: minimum number of shared barcodes\n");
         printf("       read-s1 (8)     - step 1: minimum number of reads per barcode\n");
         printf("       read-s2 (10)    - step 2: minimum number of reads per barcode\n");
         printf("       link-s1  (8)     - step 1: minimum number of shared barcodes\n");
         printf("       link-s2  (10)    - step 2: minimum number of shared barcodes\n");
         printf("       edge     (50000) - edge length to consider for scaffolding\n");
         printf("       score    (20)    - minimum average mapping score on a barcode covered area\n");
         printf("       block    (50000) - length to determine for nearest neighbours\n");
         printf("                (0)     - contigs were produced from short reads such as Illumina\n");
         printf("       file     (0)     - do not output the sam file in order to save disk space\n");
         printf("                (1)     - sam file from bwa is saved\n");
         printf("       plot     (barcode_lengtg.png) - output image file on barcode length distributions\n");
         printf("       mkdup    (Dupmarked.bam)      - output bam file with duplicated read removed \n");
         exit(1);
    }

    memset(lenplot10x,'\0',200);
    memset(mkdup,'\0',200);
    m_score = 50;
    n_nodes = 20;
    n_debug = 1;

    if(mkdup_flag == 1)
      file_tag = 1;
    strcpy(toolname,"bwa");
    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-gap"))
       {
         sscanf(argv[++i],"%d",&gap_len); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-nodes"))
       {
         sscanf(argv[++i],"%d",&n_nodes);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-sam"))
       {
         run_align = 0;
         sam_flag = 1;
         sscanf(argv[++i],"%s",samname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-bam"))
       {
         run_align = 0;
         sam_flag = 2;
         sscanf(argv[++i],"%s",bamname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-data"))
       {
         run_align = 1;
         file_tag = 2;
         sam_flag = 3;
         sscanf(argv[++i],"%s",datname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-noalign"))
       {
         run_align = 0;
         file_tag = 2;
         sam_flag = 3;
         sscanf(argv[++i],"%s",datname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-align"))
       {
         memset(toolname,'\0',500);
         sscanf(argv[++i],"%s",toolname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-mkdup"))
       {
         memset(toolname,'\0',500);
         sscanf(argv[++i],"%s",mkdup);
	 mkdup_flag = 1;
         file_tag = 1;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%d",&mscore);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-longread"))
       {
         sscanf(argv[++i],"%d",&n_longread);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-plot"))
       {
         plot_10x = 1;
         sscanf(argv[++i],"%s",lenplot10x);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-rc"))
       {
         sscanf(argv[++i],"%d",&rc_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-matrix"))
       {
         sscanf(argv[++i],"%d",&len_matrx);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-edge"))
       {
         sscanf(argv[++i],"%d",&len_edges);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-help"))
       {
         printf("Usage: %s -nodes 30 -matrix 2000 -reads 4 -score 20 -edge 50000 -link 3 -block 50000 -mkdup Dupmarked.bam -plot barcode_lengtg.png <input_assembly_fasta/q_file> <Input_read_1>> <Input_read_2> <Output_scaffold_file>\n",argv[0]);
         printf("       nodes    (30)    - number of CPUs requested\n");
         printf("       matrix   (2000)  - relation matrix size\n");
         printf("       reads    (10)    - step 1 and 2: minimum number of reads per barcode\n");
         printf("       link     (8)     - step 1 and 2: minimum number of shared barcodes\n");
         printf("       read-s1 (8)     - step 1: minimum number of reads per barcode\n");
         printf("       read-s2 (10)    - step 2: minimum number of reads per barcode\n");
         printf("       link-s1  (8)     - step 1: minimum number of shared barcodes\n");
         printf("       link-s2  (10)    - step 2: minimum number of shared barcodes\n");
         printf("       edge     (50000) - edge length to consider for scaffolding\n");
         printf("       score    (20)    - minimum average mapping score on a barcode covered area\n");
         printf("       block    (50000) - length to determine for nearest neighbours\n");
         printf("                (0)     - contigs were produced from short reads such as Illumina\n");
         printf("       file     (0)     - do not output the sam file in order to save disk space\n");
         printf("                (1)     - sam file from bwa is saved\n");
         printf("       plot     (barcode_lengtg.png) - output image file on barcode length distributions\n");
         printf("       mkdup    (Dupmarked.bam)      - output bam file with duplicated read removed \n");
         exit(1);
       }
       else if(!strcmp(argv[i],"-block"))
       {
         sscanf(argv[++i],"%d",&len_block);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&n_debug);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-reads"))
       {
         sscanf(argv[++i],"%d",&num_reads);
         num1_reads = num_reads;
         num2_reads = num_reads; 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-read-s1"))
       {
         sscanf(argv[++i],"%d",&num1_reads); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-read-s2"))
       {
         sscanf(argv[++i],"%d",&num2_reads); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-link"))
       {
         sscanf(argv[++i],"%d",&num_links);
         num1_links = num_links;
         num2_links = num_links; 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-link-s1"))
       {
         sscanf(argv[++i],"%d",&num1_links); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-link-s2"))
       {
         sscanf(argv[++i],"%d",&num2_links); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-uplink"))
       {
         sscanf(argv[++i],"%d",&uplinks); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_tag); 
         args=args+2;
       }
    }

    pid = getpid();
    memset(tempa,'\0',2000);
    if (!getcwd(tempa, sizeof(tempa)))
    {
      exit(1);
    } 
    memset(tmpdir,'\0',2000);
    memset(workdir,'\0',2000);
    sprintf(tmpdir,"%s/",tempa);
    sprintf(workdir,"%s/tmp_rununik_%d/",tempa,pid);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mkdir %s",workdir);

    RunSystemCommand(syscmd);

    if(chdir(workdir) == -1)
    {
      printf("System command error: chdir\n");
      exit(EXIT_FAILURE);
    }
     
    st = argv[0];
    ed = strrchr(argv[0],'/');
    memset(tempc,'\0',2000);
    strncpy(tempc,argv[0],ed-st);
    memset(bindir,'\0',2000);
    sprintf(bindir,"%s/scaff-bin",tempc);

    memset(file_tarseq,'\0',2000);
    memset(file_read1,'\0',2000);
    memset(file_read2,'\0',2000);
    memset(file_scaff,'\0',2000);
    memset(file_sfagp,'\0',2000);
    memset(file_datas,'\0',2000);
    memset(file_cover,'\0',2000);
    memset(file_plot10x,'\0',2000);
    memset(file_mkdup,'\0',2000);

    sprintf(file_tarseq,"%s/%s",tempa,argv[args]);
    if((run_align == 0)||(file_tag == 2))
    {
      sprintf(file_scaff,"%s/%s",tempa,argv[args+1]);
      sprintf(file_sfagp,"%s/%s.agp",tempa,argv[args+1]);
      sprintf(file_cover,"%s/%s.cov",tempa,argv[args+1]);
    }
    else
    {
      sprintf(file_read1,"%s/%s",tempa,argv[args+1]);
      sprintf(file_read2,"%s/%s",tempa,argv[args+2]);
      sprintf(file_scaff,"%s/%s",tempa,argv[args+3]);
      sprintf(file_sfagp,"%s/%s.agp",tempa,argv[args+3]);
      sprintf(file_cover,"%s/%s.cov",tempa,argv[args+3]);
    }

    if(n_longread ==0)
      mscore = 1;  
    if((namef = fopen(file_tarseq,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args]);
        exit(1);
      }
      else
      {
        memset(file_tarseq,'\0',2000);
        strcpy(file_tarseq,argv[args]);
        printf("Input target assembly file1: %s\n",file_tarseq);
      }
    }
    else
    {
      printf("Input target assembly file2: %s\n",file_tarseq);
    } 

    if(run_align)
    {
      if(file_tag != 2)
      {
        if((namef = fopen(file_read1,"r")) == NULL)
        {
          printf("File not in the working directory!\n");
          if((namef = fopen(argv[args+1],"r")) == NULL)
          {
            printf("File %s not found and please copy it to your working directory!\n",argv[args+1]);
            exit(1);
          }
          else
          {
            memset(file_read1,'\0',2000);
            strcpy(file_read1,argv[args+1]);
            printf("Input read1 file: %s\n",file_read1);
          }
        }
        else
        {
          printf("Input read1 file: %s\n",file_read1);
        } 

        if((namef = fopen(file_read2,"r")) == NULL)
        {
          printf("File not in the working directory!\n");
          if((namef = fopen(argv[args+2],"r")) == NULL)
          {
            printf("File %s not found and please copy it to your working directory!\n",argv[args+2]);
            exit(1);
          }
          else
          {
            memset(file_read2,'\0',2000);
            strcpy(file_read2,argv[args+2]);
            printf("Input read2 file: %s\n",file_read2);
          }
        }
        else
        {
          printf("Input read2 file: %s\n",file_read2);
        }
      }
      else 
      {
        sprintf(file_datas,"%s/%s",tempa,datname);
          printf("www: %s %s\n",file_datas,datname);
        if((namef = fopen(file_datas,"r")) == NULL)
        {
          printf("File not in the working directory!\n");
          if((namef = fopen(datname,"r")) == NULL)
          {
            printf("File %s not found and please copy it to your working directory!\n",datname);
            exit(1);
          }
          else
          {
            memset(file_datas,'\0',2000);
            strcpy(file_datas,datname);
            printf("Input data file: %s\n",file_datas);
          }
        }
        else
        {
          printf("Input data file: %s\n",file_datas);
        }
      }
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_screeN %s cleaN.fasta > try.out",bindir,file_tarseq);

    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_fastq -name tarseq -len 10 cleaN.fasta tarseq.fastq tarseq.tag > try.out",bindir);

    RunSystemCommand(syscmd);
    
    if((run_align)&&(strcmp(toolname,"bwa") == 0))
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/bwa index tarseq.fastq > try.out",bindir);

      RunSystemCommand(syscmd);

      if(file_tag == 0)
      {
        memset(syscmd,'\0',2000);
        if(n_longread == 0)
          sprintf(syscmd,"%s/bwa mem -t %d tarseq.fastq %s %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,n_nodes,file_read1,file_read2,"($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"%s/bwa mem -t %d tarseq.fastq %s %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,n_nodes,file_read1,file_read2,"($2<100)&&($5>=0){print $1,$2,$3,$4,$5}");
      }
      else if(file_tag == 1)
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa mem -t %d tarseq.fastq %s %s > align.sam",bindir,n_nodes,file_read1,file_read2);
    
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        if(n_longread ==0)
          sprintf(syscmd,"cat align.sam | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat","($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"cat align.sam | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat","($2<100)&&($5>=0){print $1,$2,$3,$4,$5}");
      }
      else if(file_tag == 2)
      {
        memset(syscmd,'\0',2000);
        if(n_longread ==0)
          sprintf(syscmd,"%s/scaff_FilePreProcess -t 2 -n 1 %s - |%s/bwa mem -p -t %d tarseq.fastq -  | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,file_datas,bindir,n_nodes,"($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"%s/scaff_FilePreProcess -t 2 -n 1 %s - |%s/bwa mem -p -t %d tarseq.fastq -  | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,file_datas,bindir,n_nodes,"($2<100)&&($5>=0){print $1,$2,$3,$4,$5}");
      }
      else 
      {
        printf("File input error!\n");
        exit(1);
      }  
      printf("%s\n",syscmd);
    
      RunSystemCommand(syscmd);
    }
    else if((run_align)&&(strcmp(toolname,"smalt") == 0))
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/smalt index -k 19 -s 11 hash_genome tarseq.fastq  > try.out",bindir);
    
      RunSystemCommand(syscmd);
      if(file_tag == 0)
      {
        memset(syscmd,'\0',2000);
        if(n_longread ==0)
          sprintf(syscmd,"%s/smalt map -i 1200 -j 20 -m 30 -r 888 -f samsoft -n %d -O hash_genome %s %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,n_nodes,file_read1,file_read2, "($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"%s/smalt map -i 1200 -j 20 -m 30 -r 888 -f samsoft -n %d -O hash_genome %s %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,n_nodes,file_read1,file_read2, "($2<100)&&($5>=0){print $1,$2,$3,$4,$5}");
      }
      else
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/smalt map -i 1200 -j 20 -m 30 -r 888 -f samsoft -n %d -o align.sam -O hash_genome %s %s > try.out",bindir,n_nodes,file_read1,file_read2);
    
        RunSystemCommand(syscmd);
        memset(syscmd,'\0',2000);
        if(n_longread ==0)
          sprintf(syscmd,"cat align.sam | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat","($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"cat align.sam | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat","($2<100)&&($5>=0){print $1,$2,$3,$4,$5}");
      }
      printf("%s\n",syscmd);
    
      RunSystemCommand(syscmd);
    }
    else if(run_align == 0)
    {
      if(sam_flag == 1)
      {
        memset(syscmd,'\0',2000);
        if(n_longread ==0)
          sprintf(syscmd,"cat %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",samname,"($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"cat %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",samname,"($2<100)&&($5>=0){print $1,$2,$3,$4,$5}");
        printf("%s\n",syscmd);
    
        RunSystemCommand(syscmd);
      }
      else if(sam_flag == 2)
      {
        memset(syscmd,'\0',2000);
//          sprintf(syscmd,"%s/samtools view %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,bamname,"($2<100)&&($5>0){print $1,$2,$3,$4,$5}");
        if(n_longread ==0)
          sprintf(syscmd,"%s/samtools view %s | awk '%s' | egrep -v '^@' > align0.dat",bindir,bamname,"($4!=0)&&($2<100)&&($5>0){print $1,$12,$2,$3,$4,$5}");
        else
          sprintf(syscmd,"%s/samtools view %s | awk '%s' | egrep -v '^@' > align0.dat",bindir,bamname,"($4!=0)&&($2<100)&&($5>=0){print $1,$12,$2,$3,$4,$5}");
        printf("%s\n",syscmd);
    
        RunSystemCommand(syscmd);
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/scaff_bwa-barcode tarseq.tag align0.dat align.dat > try.out",bindir);
        printf("%s\n",syscmd);
        RunSystemCommand(syscmd);
      }
      else if(sam_flag == 3)
      {
        memset(syscmd,'\0',2000);
//        sprintf(syscmd,"cp %s align.dat",datname);
        sprintf(syscmd,"cp %s align.dat",datname);
        printf("%s\n",syscmd);
        RunSystemCommand(syscmd);
      }
      else
      {
        printf("File input error:\n");
        exit(1);
      }

    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cp %s/BC.code BC.code",bindir);
    RunSystemCommand(syscmd);
      
    memset(syscmd,'\0',2000);
    printf("%s/scaff_bwa -edge %d tarseq.tag BC.code align.dat align2.dat > try.out",bindir,len_edges);
    sprintf(syscmd,"%s/scaff_bwa -edge %d -mod %d tarseq.tag BC.code align.dat align2.dat > try.out",bindir,len_edges,rc_flag);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_barcode-sort align2.dat align.sort > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_contigs-sort align.sort align.sort2 > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_mapping-sort -block %d -reads %d -score %d align.sort2 barcodes.clust > try.out",bindir,len_block,num1_reads,mscore);

    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_mapping-clean barcodes.clust barcodes.clean > try.out",bindir);
    
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_matrix -file 1 -matrix %d -link %d -uplink %d -longread %d barcodes.clean tarseq.tag contig.dat > scaff.out",bindir,len_matrx,num1_links,uplinks,n_longread);
    
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_output -longread %d -gap %d tarseq.fastq contig.dat genome.fastq genome.agp > try.out",bindir,n_longread,gap_len);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_RDplace align.sort2 genome.agp align.sort3 > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_mapping-sort -block %d -reads %d -score %d align.sort3 barcodes.clust2 > try.out",bindir,len_block,num2_reads,mscore);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_mapping-clean barcodes.clust2 barcodes.clean2 > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_fastq -name scaff10x -len 1 genome.fastq genome-new.fastq genome.tag > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_matrix -file 1 -matrix %d -link %d -uplink %d -longread %d barcodes.clean2 genome.tag contig.dat2 > scaff.out2",bindir,len_matrx,num2_links,uplinks,n_longread);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_output -longread %d -gap %d genome.fastq contig.dat2 genome2.fastq genome2.agp > try.out",bindir,n_longread,gap_len);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_rename genome2.fastq genome.fasta > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mv genome.fasta %s",file_scaff);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scaff_agp2agp genome.agp genome2.agp tarseq.tag genome-all.agp > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mv genome-all.agp %s",file_sfagp);
    RunSystemCommand(syscmd);

    if(mkdup_flag == 1)
    {
      sprintf(file_mkdup,"%s/%s",tempa,mkdup);
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools view -bS align.sam > Sorted_names.bam ",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools fixmate -m Sorted_names.bam Fixmate.bam > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools sort -o Sorted.bam Fixmate.bam > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf Sorted_names.bam Fixmate.bam ");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/samtools markdup -r -s Sorted.bam Dupmarked.bam > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv Dupmarked.bam %s",file_plot10x);
      RunSystemCommand(syscmd);
    }
    if(plot_10x == 1)
    {
      sprintf(file_plot10x,"%s/%s",tempa,lenplot10x);
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_bwa -edge 300000000 -mod %d tarseq.tag BC.code align.dat align.size > try.out",bindir,rc_flag);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep ^A align.size > align2.dat_AAA");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep ^C align.size > align2.dat_CCC");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep ^G align.size > align2.dat_GGG");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep ^T align.size > align2.dat_TTT");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.size");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_barcode-sort align2.dat_AAA align.sort_AAA > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align2.dat_AAA");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_barcode-sort align2.dat_CCC align.sort_CCC > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align2.dat_CCC");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_barcode-sort align2.dat_GGG align.sort_GGG > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align2.dat_GGG");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_barcode-sort align2.dat_TTT align.sort_TTT > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align2.dat_TTT");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort align.sort_AAA align.sort2_AAA > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort_AAA");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort align.sort_CCC align.sort2_CCC > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort_CCC");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort align.sort_GGG align.sort2_GGG > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort_GGG");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort align.sort_TTT align.sort2_TTT > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort_TTT");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort -break 1 align.sort2_AAA align.sort3_AAA > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort2_AAA");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort -break 1 align.sort2_CCC align.sort3_CCC > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort2_CCC");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort -break 1 align.sort2_GGG align.sort3_GGG > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort2_GGG");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_contigs-sort -break 1 align.sort2_TTT align.sort3_TTT > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort2_TTT");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cat align.sort3_AAA align.sort3_CCC align.sort3_GGG align.sort3_TTT > align.size4");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_PCRdup align.size4 align.size5 > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.sort3_AAA align.sort3_CCC align.sort3_GGG align.sort3_TTT");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.size4 ");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_barcode-length -reads 3 align.size5 align.length-5 > try.out",bindir);
      RunSystemCommand(syscmd);

/*
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf align.size5 ",bindir);
      RunSystemCommand(syscmd);
                                        */
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cat align.length-5 | awk '{print $1,$3}' | sort -n -k 2 > sample-bcl-5.dat");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_lengthdis sample-bcl-5.dat sample-10x.freq > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/human-bcl-5.freq human-bcl-5.freq",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/hummingbird-bcl-5.freq hummingbird-bcl-5.freq",bindir);
      RunSystemCommand(syscmd);
      
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/fAnaTes1-bcl-5.freq fAnaTes1-bcl-5.freq",bindir);
      RunSystemCommand(syscmd);
      
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/fSimDai1-bcl-5.freq fSimDai1-bcl-5.freq",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"bash %s/plot-10x-length.sh",bindir);
      RunSystemCommand(syscmd);
      
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp plot10x.png %s",lenplot10x);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv %s %s",lenplot10x,file_plot10x);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"sort -k 2,2 -k 5,5n align.length-5 > align.length-sort");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/scaff_barcode-cover align.length-sort break.dat cover.dat > break.out",bindir);
      RunSystemCommand(syscmd);
    
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp cover.dat %s",file_cover);
      RunSystemCommand(syscmd);
    
    }    
    if(n_debug == 0)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf * > /dev/null");
    
      RunSystemCommand(syscmd);
      chdir(tmpdir);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf %s > /dev/null",workdir);
    
      RunSystemCommand(syscmd);
    }
    return EXIT_SUCCESS;

}
/* end of the main */



/*   subroutine to sort out read pairs    */
/* =============================== */
void Read_Index(int nSeq, char *namefile)
/* =============================== */
{
     int i,j,nseq;
     int i_reads,n_reads,c_reads,insertsize=0;
     FILE *namef;
     char *ptr;
     char line[500];
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);


     if((namef = fopen(namefile,"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     nseq = 0;
     i_reads = 0;
     n_reads = 0;
     c_reads = 0;
     while(!feof(namef))
     {
       if(fgets(line,500,namef) == NULL)
            printf("Data input file problem!\n");
       if(feof(namef)) break;
       nseq++;
     }
     fclose(namef);

     nseq = 2*nseq;
     if((ctg_head = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_head\n");
       exit(1);
     }
     if((ctg_list = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((core_list = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - core_list\n");
       exit(1);
     }
     if((read2contig = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - read2contig\n");
       exit(1);
     }

     nseq = nseq*3;
     S_Name=cmatrix(0,nseq,0,Max_N_NameBase);
     if((insert_siz = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - insert\n");
       exit(1);
     }
     if((insert_dev = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - insert\n");
       exit(1);
     }
 
     if((namef = fopen(namefile,"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

     j = 0;
     insertsize = 0;
     while(!feof(namef))
     {
       int nPair=0,len=0;
       char line2[500],line3[500],base[500];

       if(fgets(line,500,namef) == NULL)
            printf("Data input file problem!\n");
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       insertsize = 0;
       if((strncmp(line,"readnames",9))==0)
       {
         i=0;
         c_reads = 0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
                memset(base,'\0',500);
//                len=strlen(ptr);
//                strncpy(base,ptr,len-1);
                strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
//       printf("creads: %d %d\n",c_reads,n_reads);
         if(n_group>0)
           ctg_list[n_group-1]=n_reads;
         n_group++;
         n_reads = 0;
       }
       else
       {      
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line3," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(nPair>1)
            {
              if(i==(nPair-2))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                insertsize = atoi(base);
              }
            }
         }
         i=0;
         j=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {

              len=strlen(ptr);
              if(nPair==1)
                strncpy(S_Name[i_reads],ptr,len-1);
              else
                strncpy(S_Name[i_reads],ptr,len);
//       printf("reads: %d %d %d %d %s\n",j,i_reads,insertsize,c_reads,S_Name[i_reads]);
              i_reads++;
              j++;
            }
//            else if(insertsize<50000)
            else if((insertsize<50000)&&(c_reads>4))
            {
              if(i==1)
              {
                len=strlen(ptr);
                strncpy(S_Name[i_reads],ptr,len);
                i_reads++;
                j++;
              }
              else if((i==2)&&(i<(nPair-2)))
              {
                len=strlen(ptr);
                strncpy(S_Name[i_reads],ptr,len);
                i_reads++;
                j++;
              }
/*              else if(i==(nPair-2))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                for(k=0;k<j;k++)
                   insert_siz[i_reads-k-1] = atoi(base);
              }
              else if(i==(nPair-1))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                for(k=0;k<j;k++)
                   insert_dev[i_reads-k-1] = atoi(base);
              }   */
            }
         }
         n_reads = n_reads+j;
       }
       c_reads++;
     }
     fclose(namef);
     printf("contig: %d %d\n",n_reads,i_reads);
     ctg_list[n_group-1]=n_reads;
     ctg_head[0]=0;
     for(i=1;i<n_group;i++)
        ctg_head[i]=ctg_head[i-1]+ctg_list[i-1];
      
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   to swap the string arrays           */
/* ============================================= */
void s_swap2(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase2];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String2(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase2];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap2(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap2(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap2(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

