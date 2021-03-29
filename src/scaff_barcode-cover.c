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
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 50000 
#define Max_N_NameBase 50
static char **S_Name;
static int *hit_locus1,*hit_locus2,*hit_length,*hit_barlen,*hit_mscore,*hit_reads;

/* SSAS default parameters   */
static int IMOD=0;
static int n_type=0;
static int barreads=10;
static int file_flag=2;
static int break_flag=0;
static int mp_score=20;
static int nContig=0;
static int max_len = 0;
static int read_len = 150;
static int num_contig = 0;
static int min_ratio = 15;
static int min_cover = 50;

fasta *expt;

static char rc_char[500000];
static char rc_sub[5000];

int ReverseComplement(int seqdex)
{
        int i,len;
        char *tp,*dp;
        fasta *seqp;

        seqp=expt+seqdex;
        len=seqp->length;
        memset(rc_sub,'\0',5000);
        dp=rc_sub;      
        tp = seqp->data+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int Reverse_Complement_Contig(char c_array[],int num_len)
{
        int i,len;
        char *tp,*dp;

        len=num_len;
        dp=rc_char;
        tp = c_array+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int main(int argc, char **argv)
{
    FILE *namef;
    int i,j,nSeq,args,idt;
    int n_contig,n_reads,n_readsMaxctg,nseq;
    fasta *seq;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Readname_match(fasta *seq,char **argv,int args,int nSeq,int nRead);
    void Barcode_Process(char **argv,int args,int nSeq);
    void Memory_Allocate(int arr);
    char line[2000]={0},tempc1[60],cc[60],RC[2],readname[60],*st,*ed;
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    void Read_Pairs(char **argv,int args,fasta *seq,int nSeq);

    seq=NULL;
    if(argc < 2)
    {
      printf("Usage: %s -score 20 -cover 50 -ratio 15 <input_barcode_file> <barcode_sorted file>\n",argv[0]);

      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&read_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%d",&mp_score);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-ratio"))
       {
         sscanf(argv[++i],"%d",&min_ratio);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cover"))
       {
         sscanf(argv[++i],"%d",&min_cover);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-break"))
       {
         sscanf(argv[++i],"%d",&break_flag);
         args=args+2;
       }
    }

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      if(fgets(line,2000,namef) == NULL)
      {
//        printf("fgets command error:\n);
      }	    
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 
   
/*
    nRead=0;
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: args+1 \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef);   */ 

    if((hit_locus1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus\n");
      exit(1);
    }
    if((hit_locus2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus\n");
      exit(1);
    }
    if((hit_length = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus2\n");
      exit(1);
    }
    if((hit_barlen = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_barlen\n");
      exit(1);
    }
    if((hit_mscore = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus2\n");
      exit(1);
    }
    if((hit_reads = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_reads\n");
      exit(1);
    }

    nSeq=nseq;
    S_Name=cmatrix(0,nseq,0,Max_N_NameBase);

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    num_contig = 0;
    while(fscanf(namef,"%s %s %d %d %d %d %d %s %d",cc,S_Name[i],&hit_barlen[i],&hit_length[i],&hit_locus1[i],&hit_locus2[i],&hit_reads[i],cc,&hit_mscore[i])!=EOF)
    {
        ed = strrchr(S_Name[i],'_');
        idt = atoi(ed+1);
        if(idt >= num_contig)
          num_contig = idt; 
        if(hit_length[i] > max_len)
          max_len = hit_length[i];
        i++;
    }
    fclose(namef);


    n_reads=i;
    Barcode_Process(argv,args,n_reads);

    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Barcode_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n,g_size,max_barlen;
     int num_hits,*hit_bccover,*hit_rdcover,*hit_ctglens,*cov_genome;
     FILE *namef;
     float ave_barlen,R50,R60,R70,R80,R90;
     long num_cover,ave_cover,t_SQbases,t_BCbases,t_RDbases,totalBarlen,totalHalf;
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     void ArraySort_Int2(int n, int *arr, int *brr);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     int N50,N60,N70,N80,N90,M50,M60,M70,M80,M90;


     g_size = max_len+1000;;      
     if((cov_genome= (int *)calloc(g_size,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - hit_maps\n");
       exit(1);
     }

     num_contig = num_contig +1;          
     if((hit_bccover = (int *)calloc(num_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_bccover\n");
       exit(1);
     }
     if((hit_rdcover = (int *)calloc(num_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_rdcover\n");
       exit(1);
     }
     if((hit_ctglens = (int *)calloc(num_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - hit_ctglens\n");
       exit(1);
     }

     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     num_hits =0;

     t_SQbases = 0;
     t_BCbases = 0;
     t_RDbases = 0;
 
     for(i=0;i<(nSeq-1);i++)
     {
        int stopflag=0;
        int idt;
        char *ed;

        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(S_Name[i],S_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        ed = strrchr(S_Name[i],'_');
        idt = atoi(ed+1);
        hit_ctglens[idt] = hit_length[i];
        num_hits = j-i;
        if((num_hits>=5)&&(hit_length[i]> 50000)) 
        {
          int ctg_len = hit_length[i];
          int set_cover = 0;
          int n_reads = 0;
          int n_bases = 0;
          int read_cover = 0;
          int ii,jj,stopflag,mini_hit,mini_los;

          stopflag = 0;
          for(k=0;k<ctg_len;k++)
             cov_genome[k] = 0;
//          memset(cov_genome,0,4*ctg_len);
//          printf("Scaffold: 2222 %s %d %d %d\n",S_Name[i],i,num_hits,hit_length[i]);
          num_cover = 0;
          ave_cover = 0;
          for(n=i;n<j;n++)
          {
             if(hit_mscore[n] >= mp_score)
             {
               for(ii=hit_locus1[n];ii<=hit_locus2[n];ii++)
               {
                  cov_genome[ii]++;
                  num_cover++;
               }
             }
             n_reads = n_reads + hit_reads[n];
//             fprintf(namef,"%s %s %d %s %d\n",S_Name[n],R_Name[idt],hit_locus[idt],M_Name[idt],hit_length[idt]);
          }
          ave_cover = num_cover;
          ave_cover = ave_cover/ctg_len;
          set_cover = 0.3*ave_cover;
          n_bases = read_len+read_len;
          n_bases = n_bases*n_reads;
          read_cover = n_bases/ctg_len;
          t_SQbases = t_SQbases + ctg_len;
          t_BCbases = t_BCbases + num_cover;
          t_RDbases = t_RDbases + n_bases;
          hit_bccover[idt] = ave_cover;
          hit_rdcover[idt] = read_cover;
          printf("Scaffold: %s %d %ld %d %d %d %d\n",S_Name[i],ctg_len,ave_cover,set_cover,n_reads,n_bases,read_cover);
/*          if(strcmp(S_Name[i],"tarseq_205")==0)
          {
            for(ii=0;ii<ctg_len;ii++)
                 printf("xxx: %d %d %ld %d\n",ii,cov_genome[ii],ave_cover,set_cover);
          }  */
          for(n=0;n<ctg_len;n++)
          {
             if(n<50000)
               cov_genome[n] = -1;
             else if(n>(ctg_len-50000))
               cov_genome[n] = -1;
             else if(cov_genome[n] >= min_cover)
             {
               if(cov_genome[n] >= set_cover)
                 cov_genome[n] = -1; 
             } 
          }

          for(ii=0;ii<ctg_len;ii++)
          {
             jj = ii+1;
             stopflag = 0;
             while((jj<ctg_len)&&(stopflag==0))
             {
               if(cov_genome[jj]>=0)
               {
                 jj++;
               }
               else
                 stopflag=1;
             }
             if((jj-ii) > 20)
             {
               mini_hit = 99999;
               mini_los = 0;
               for(n=(ii+1);n<(jj-1);n++)
               {
                  if(cov_genome[n] < mini_hit)
                  {
                    mini_hit = cov_genome[n];
                    mini_los = n;
                  }
               }
//               if(mini_hit < min_cover)
               set_cover = min_ratio*ave_cover;
               set_cover = set_cover/100;
//                 printf("www: %d %d %d %ld %d\n",mini_los,ctg_len,mini_hit,ave_cover,set_cover);
               if(mini_hit < set_cover)
               {
                 int idt;
                 char *ed;
                 ed = strrchr(S_Name[i],'_');
                 idt = atoi(ed+1); 
                 printf("Break: %d %d %d %d %ld %d\n",idt,mini_los,ctg_len,mini_hit,ave_cover,set_cover);
                 fprintf(namef,"Break: %d %d %d %d %ld\n",idt,mini_los,ctg_len,mini_hit,ave_cover);
               }
             }
             ii = jj-1;
          }
        }
//        else
//          fprintf(namef,"%s %s %d F %d\n",S_Name[i],R_Name[i],hit_locus[i],hit_length[i]);         
        i = j - 1; 
     }
     fclose(namef);

     if((namef = fopen(argv[args+2],"w")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

     totalBarlen = 0;
     for(i=0;i<nSeq;i++)
     {
        hit_reads[i] = i;
        totalBarlen = totalBarlen + hit_barlen[i];
     }

     ArraySort2_Int2(nSeq,hit_barlen,hit_reads);
     max_barlen = hit_barlen[0];
     ave_barlen = totalBarlen/nSeq;
     totalHalf = 0;
     N50 = 0;
     N60 = 0;
     N70 = 0;
     N80 = 0;
     N90 = 0;
     fprintf(namef,"Barcode length: n = %d, ave = %f, largest = %d\n",nSeq,ave_barlen,max_barlen);

     for(i=0;i<nSeq;i++)
     {
        totalHalf = totalHalf+hit_barlen[i];
        if((totalHalf >= 0.5*totalBarlen)&&(N50==0))
        {
          M50 = 0.5*totalBarlen/(i+1);
          R50 = M50;
          R50 = R50/hit_barlen[i];
          fprintf(namef,"        N50 = %d, n = %d\n",hit_barlen[i],i+1);
          N50 = 1;
        }
        if((totalHalf >= 0.6*totalBarlen)&&(N60==0))
        {
          M60 = 0.6*totalBarlen/(i+1);
          R60 = M60;
          R60 = hit_barlen[i]/R60;
          R60 = R60*R50;
          fprintf(namef,"        N60 = %d, n = %d\n",hit_barlen[i],i+1);
          N60 = 1;
        }
        if((totalHalf >= 0.7*totalBarlen)&&(N70==0))
        {
          M70 = 0.7*totalBarlen/(i+1);
          R70 = M70;
          R70 = hit_barlen[i]/R70;
          R70 = R70*R50;
          fprintf(namef,"        N70 = %d, n = %d\n",hit_barlen[i],i+1);
          N70 = 1;
        }
        if((totalHalf >= 0.8*totalBarlen)&&(N80==0))
        {
          M80 = 0.8*totalBarlen/(i+1);
          R80 = M80;
          R80 = hit_barlen[i]/R80;
          R80 = R80*R50;
          fprintf(namef,"        N80 = %d, n = %d\n",hit_barlen[i],i+1);
          N80 = 1;
        }
        if((totalHalf >= 0.9*totalBarlen)&&(N90==0))
        {
          M90 = 0.9*totalBarlen/(i+1);
          R90 = M90;
          R90 = hit_barlen[i]/R90;
          R90 = R90*R50;
          fprintf(namef,"        N90 = %d, n = %d\n",hit_barlen[i],i+1);
          N90 = 1;
          i = nSeq;
        }
     }

     fprintf(namef,"        N100 = %d, n = %d\n",hit_barlen[nSeq-1],nSeq);

     fprintf(namef,"Barcode cover: %ld\n",t_BCbases/t_SQbases);
     fprintf(namef,"Read cover: %ld\n",t_RDbases/t_SQbases);

     for(i=0;i<num_contig;i++)
     {
        fprintf(namef,"Cover: tarseq_%d %d %d %d\n",i,hit_ctglens[i],hit_bccover[i],hit_rdcover[i]);
     }
     fclose(namef);

}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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

