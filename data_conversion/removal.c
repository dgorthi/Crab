#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include "crab_databuf.h"
#include "filterbank.h"
#include "rand_num_generator.h"
#include <dirent.h>

#define TRUE 1
#define FALSE 0 
#define REMOVE_BAD_CHANNEL TRUE
#define REPLACE_BAD_ZERO FALSE
#define REMOVE_NARROWBAND_RFI TRUE
#define REPLACE_NARROW_ZERO FALSE
#define REMOVE_BROADBAND_RFI FALSE
#define REPLACE_BROAD_ZERO FALSE
#define HIGHPASS_FILTER FALSE

void remove_bad_channel(short int data[], int data_len, int start_channel, int stop_channel);

int calc_stats_broad(double data[], short int replace[], double avg[], double sd[], int num_samples, int start_channel_offset, int start_channel, int stop_channel, int window_size, int channels_per_spec, long start_time);

int calc_stats_narrow(short int data[], short int replace[], double avg[], double sd[], int data_len, int start_channel_offset, int start_channel, int stop_channel, int window_size, int channels_per_spec, long start_time);

static int myCompare (const void * a, const void * b);

void sort(const char *arr[], int n);

int main(){
	
	char f_fil[256];
	int start_channel;
	int stop_channel;
	double sum_channels;
	int n_data_points;
	long seed;
	int counter;
	int num_channel_per_spec = 4096;
	long long chunk;
	char temp[255];
	int num_files = 0;
	struct tm        *now;
	time_t           rawtime;
	FILE * fp;
	FILE * crab_file;
	char fileName[255];
	sum_channels = 0;
	seed = 5;
/*
	DIR *p;
	struct dirent *pp;     
	p = opendir ("./");
	if (p != NULL)
	{
		while ((pp = readdir (p))!=NULL) {
		int length = strlen(pp->d_name);
	      		if (strncmp(pp->d_name, "data", 4) == 0 && strncmp(pp->d_name+length-4, ".fil", 4) != 0) {
				num_files++;
			}
		}
	}
	rewinddir(p);
	printf("%d",num_files);
	char files[162][255];
	
	char **files;
	files = (char **)malloc(num_files * sizeof(char*));
	for (int i = 0; i < num_files; i++)
    		files[i] = (char *)malloc((255) * sizeof(char)); 
	if (!files){
		printf("malloc failed");
	}

// Clear-up
	if (p != NULL)
	{
		while ((pp = readdir (p))!=NULL) {
		int length = strlen(pp->d_name);
	      		if (strncmp(pp->d_name, "data", 4) == 0 && strncmp(pp->d_name+length-4, ".fil", 4) != 0) {
	//	  		puts (pp->d_name);
				strncpy(fileName, pp->d_name, 254);
        			fileName[254] = '\0';
	//			printf("%s\n",fileName);
				strcpy(files[counter], fileName);
	//			files[counter] = fileName;
				for(int i=0;i<= num_files;i++)
      				for(int j=i+1;j<= num_files;j++){
					 if(strcmp(files[i],files[j])>0){
					    strcpy(temp,files[i]);
					    strcpy(files[i],files[j]);
					    strcpy(files[j],temp);
					 }
				}
	//			printf("%s\n", files[counter]);
				counter++;
	      		}
	   	}
	    	(void) closedir (p);
	}
	//sort(files, 90);
	printf("%s\n", files[80]);

	for(int i = 0; i<90;i++){
		printf("%d\n", i);
		printf("%s\n\n",files[i]);
	}
*/	
	printf("\n\nopen new filterbank file...\n\n");
	time(&rawtime);
	now = localtime(&rawtime);
	strftime(f_fil,sizeof(f_fil), "./data_%Y-%m-%d_%H-%M-%S.fil",now);
	WriteHeader(f_fil);
	printf("write header done!\n");
	crab_file=fopen(f_fil,"a+");

	//for(int i = 81; i <(81+ 1); i++){
	//fp = fopen("data_2017-08-09_11-43-07.fil", "rb");
	//fp = fopen("data_2017-07-19_15-05-13.fil", "rb");
	fp = fopen("test.fil", "rb");
	//n_data = 2097152*9  # Number of data points per buffer
	//n_flag = 1024*9      # Number of flags per buffer


//	fp = fopen(files[i], "rb");
	
	fseek(fp, 0L, SEEK_END);
	chunk = (ftell(fp)-363)/2;
	//long fsize = ftell(fp);
	//int no_cpy_buf = (fsize/(2097152*18));
	//printf("no copy buf %d\n", fsize);
	//n_data_points = (fsize-363)/sizeof(short int);
	n_data_points = 246345728;  //num data points per minute of data
	printf("num data points: %d\n", n_data_points);
	fseek(fp, 363, SEEK_SET);
	//fseek (fp, 0, SEEK_SET);
	printf("starting write data...\n");	
	for(long long i = 0; i <= chunk-n_data_points; i+= n_data_points){
		int num_samples = n_data_points/num_channel_per_spec;
		
		short int *I = (short int*)malloc(sizeof(short int)*n_data_points);
		double *time_sum = (double*)malloc(sizeof(double)*num_samples);
		double *channel_avg = (double*)malloc(sizeof(double)*n_data_points);
		double *channel_SD = (double*)malloc(sizeof(double)*n_data_points);
		double *channel_avg_total = (double*)malloc(sizeof(double)*num_samples);
		double *channel_sd_total = (double*)malloc(sizeof(double)*num_samples);
		if (!I || !channel_avg || !channel_SD || !channel_avg_total || !channel_sd_total){
			printf("Malloc Failed");
		}
		int size = fread(I, sizeof(short int), n_data_points, fp);
		start_channel = 1800;
		stop_channel = 2730;
		//start_channel = 240;
		//stop_channel = 1264;

		
		for(int j = 0; j < n_data_points/num_channel_per_spec; j++){
			for(int k = 0; k< num_channel_per_spec; k++){
				if (k < start_channel || k >= stop_channel){
					I[num_channel_per_spec*j+k] = 0;
				}
			}
		}

		#if REMOVE_BAD_CHANNEL == TRUE
			remove_bad_channel(I, n_data_points, start_channel, stop_channel);
		#endif

		#if REMOVE_NARROWBAND_RFI == TRUE
			for(int n = start_channel; n < stop_channel; n++){
				counter += calc_stats_narrow(I, I, channel_avg, channel_SD, num_samples, 0, n, n+1, 1001, num_channel_per_spec,0);
			}
		#endif
		printf("Number of narrowband samples removed: %d\n", counter);
		
		counter = 0;
		for(int j = 0; j < num_samples; j++){
			for(int k = 0; k< num_channel_per_spec; k++){
				//printf("DATA VALUE: %d\n", I[num_channel_per_spec*j+k]);
				//printf("AVG VALUE %d\n\n", (short int)channel_avg[num_channel_per_spec*j+k]);
				sum_channels += I[num_channel_per_spec*j+k];
			}
			time_sum[j] = sum_channels;
//			printf("TIME_SUM: %f\n\n", time_sum[j]);
			sum_channels = 0;
		}
/*	
		for(int j = 0; j < n_data_points/num_channel; j++){
				printf("Channel  value at time sample %d: %.6f\n",  j, time_sum[j]);
				printf("Channel avg at time sample %d: %.6f\n",  j, channel_avg_total[j]);
				printf("Channel SD at time sample %d: %.6f\n\n\n",  j, channel_sd_total[j]);
		}
	
	
		for(int j = 5581; j < 5583; j++){
			for(int k = start_channel; k<stop_channel; k++){
				printf("Channel %d avg at time sample %d: %.6f\n", k, j, channel_avg[j*num_channel_per_spec+k]);
				printf("Channel %d SD at time sample %d: %.6f\n", k, j, channel_SD[j*num_channel_per_spec+k]);
				printf("Channel %d value at time sample %d: %d\n\n\n", k, j, I[j*num_channel_per_spec+k]);
			}
		}
*/		
		#if REMOVE_BROADBAND_RFI == TRUE
			counter = calc_stats_broad(time_sum, I, channel_avg, channel_SD, num_samples, start_channel, start_channel, stop_channel, 1001, 1, 0);
			counter += calc_stats_broad(time_sum, I, channel_avg, channel_SD, num_samples, start_channel, start_channel, stop_channel, 60143, 1, 0);
		#endif

		printf("Number of broadband samples removed: %d\n",counter);

		#if HIGHPASS_FILTER == TRUE
			for(int j = 0; j < num_samples; j++){
                        	for(int k = 0; k< num_channel_per_spec; k++){
                                //printf("DATA VALUE: %d\n", I[num_channel_per_spec*j+k]);
                                //printf("AVG VALUE %d\n\n", (short int)channel_avg[num_channel_per_spec*j+k]);
                                I[num_channel_per_spec*j+k]-= (short int)channel_avg[num_channel_per_spec*j+k];
				}
                        }
		#endif

		fwrite(I,sizeof(short int),n_data_points,crab_file);
		free(I);
		free(channel_avg);
		free(channel_SD);
		free(channel_avg_total);	
		free(channel_sd_total);
		free(time_sum);
		fclose(fp);
	}
	fclose(crab_file);
	//free(files);

	return(0);
}

void remove_bad_channel(short int data[], int data_len, int start_channel, int stop_channel){
	int num_channel = 4096;
	double sd_new;
	double sd_old;
	double m_newM;
	double m_oldM;
	double m_newS;
	double m_oldS;
	double thresh_top;
	double thresh_bottom;
	double random;
	long seed = 5;
	int num_samples = data_len/num_channel;
	
	for(int j = 0; j < num_samples; j++){
		m_oldM = data[(j)*num_channel+start_channel];
		m_oldS = 0;
		for(int n =start_channel+1; n < stop_channel; n++){
			m_newM = m_oldM + (data[j*num_channel+n]-m_oldM)/(n);
			m_newS = m_oldS + ((double)data[j*num_channel+n]-m_oldM)*((double)data[j*num_channel+n]-m_newM);
			m_oldM = m_newM;
			m_oldS = m_newS;
		}
		sd_new = sqrt(m_newS/((stop_channel-start_channel)-1));
		sd_old = sd_new;
		thresh_top = m_newM + 3*sd_new;
		thresh_bottom = m_newM-3*sd_new;
		#if REPLACE_BAD_ZERO==TRUE
			for(int n = start_channel; n < stop_channel; n++){
				if(data[j*num_channel+n] > thresh_top || data[j*num_channel+n] < thresh_bottom){
					data[j*num_channel+n] = 0;
				}	
			}
		#else
			for(int n = start_channel; n < stop_channel; n++){
				if(data[j*num_channel+n] > thresh_top || data[j*num_channel+n] < thresh_bottom){
	//				printf("1\n");
					random = gasdev(&seed);
					while(random>1 || random<-1){
						random = gasdev(&seed);
					}
					data[j*num_channel+n] = (short int)(m_newM+sd_new*random);
				}	
			}
		#endif	
	}
}

int calc_stats_broad(double data[], short int replace[], double avg[], double sd[], int num_samples, int start_channel_offset, int start_channel, int stop_channel, int window_size, int channels_per_spec, long start_time){
	double m_oldM;
	double m_newM;
	double m_oldS;
	double m_newS;
	double data_old;
	double data_new;
	double sd_new;
	double sd_old;
	double threshold_top;
	double threshold_bottom;
	double random;
	long seed = 5; 
	int counter = 0;
	long j = start_time;
	num_samples = start_time+num_samples;
	//for(long j = start_time; j < (start_time+num_samples); j++){
	while(j < num_samples){
		if (j==start_time){
			m_oldM = data[j*channels_per_spec+start_channel-start_channel_offset];
			m_oldS = 0;
			for(int k = 1; k < window_size; k++){
				m_newM = m_oldM+ ((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset]-m_oldM)/(k+1);
				m_newS = m_oldS+((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset]-m_oldM)*((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset]-m_newM);
				m_oldM = m_newM;
				m_oldS = m_newS;
			}
			sd_new = sqrt(m_newS/(window_size-1));
			sd_old = sd_new;
//			printf("sd %.6f \n", sd_new);
//			printf("mean %.6f \n", m_newM);
//			printf("data %d \n", data[j*channels_per_spec+n]);
		}
		if (j>(window_size+start_time) && j < (start_time+(num_samples-window_size))){
			data_new = data[(j-1)*channels_per_spec+start_channel-start_channel_offset];
			data_old = data[(j-window_size-1)*channels_per_spec+start_channel-start_channel_offset];
			m_newM = m_oldM - (double)data_old/window_size + (double)data_new/window_size;
				sd_new = sqrt(sd_old*sd_old + (double)(data_new - data_old)*(data_new-m_newM + data_old - m_oldM)/(window_size-1));
				m_oldM = m_newM;
				sd_old = sd_new;
		}
		//printf("sd %.6f \n", sd_new);
		//printf("mean %.6f \n", m_newM);
		//printf("number %f\n\n", data[j*channels_per_spec+start_channel-start_channel_offset]);
		//	avg[(j-start_time)*channels_per_spec+start_channel] = m_newM;
		//	sd[(j-start_time)*channels_per_spec+start_channel] = sd_new;
		threshold_top = m_newM+sd_new*3;
		threshold_bottom = m_newM-sd_new*3;
		if (data[j*channels_per_spec+start_channel-start_channel_offset] > threshold_top || data[j*channels_per_spec+start_channel-start_channel_offset] < threshold_bottom){// && j>(window_size+start_time)){
			channels_per_spec = 1536;
			counter++;
			#if REPLACE_BROAD_ZERO==TRUE
				for(int n = start_channel; n < stop_channel; n++){
					replace[j*channels_per_spec+n] = 0;
				}
			#else
				for(int n = start_channel; n < stop_channel; n++){
	//				printf("2\n");
					random = gasdev(&seed);
					while(random>1 || random<-1){
						random = gasdev(&seed);
					}
					//printf("%d  %f\n", j, sd[(j-start_time)*channels_per_spec+n]);
//					printf("%f\n", random);
					replace[j*channels_per_spec+n] = (short int) (avg[(j-start_time)*channels_per_spec+n]+sd[(j-start_time)*channels_per_spec+n]*random);
				}
			#endif	
			channels_per_spec = 1;
			}
		j++;
	}
	return counter;
}

int calc_stats_narrow(short int data[], short int replace[], double avg[], double sd[], int num_samples, int start_channel_offset, int start_channel, int stop_channel, int window_size, int channels_per_spec, long start_time){
	double m_oldM;
	double m_newM;
	double m_oldS;
	double m_newS;
	double data_new;
	double data_old;
	double sd_new;
	double sd_old;
	double threshold_top;
	double threshold_bottom;
	double random;
	//double stat_array[1001]; 
	long seed = 5; 
	int counter = 0;
	long j = start_time;
	num_samples = start_time+num_samples;
	//for(long j = start_time; j < (start_time+num_samples); j++){
	while(j < num_samples){
		if (j==start_time){
//			stat_array = data[j*channels_per_spec+start_channel-start_channel_offset];
			m_oldM = data[j*channels_per_spec+start_channel-start_channel_offset];
			m_oldS = 0;
			for(int k = 1; k < window_size; k++){
//				stat_array[k] = ((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset];
				m_newM = m_oldM+ ((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset]-m_oldM)/(k+1);
				m_newS = m_oldS+((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset]-m_oldM)*((double)data[(j+k)*channels_per_spec+start_channel-start_channel_offset]-m_newM);
				m_oldM = m_newM;
				m_oldS = m_newS;
			}
			sd_new = sqrt(m_newS/(window_size-1));
			sd_old = sd_new;
			//printf("sd %.6f \n", sd_new);
			//printf("mean %.6f \n", m_newM);
			//printf("data %d \n", data[j*channels_per_spec+start_channel]);
		}
		if (j>(window_size+start_time) && j < (start_time+(num_samples-window_size))){
			data_new = data[(j-1)*channels_per_spec+start_channel-start_channel_offset];
			data_old = data[(j-window_size-1)*channels_per_spec+start_channel-start_channel_offset];
			m_newM = m_oldM - (double)data_old/window_size + (double)data_new/window_size;
				sd_new = sqrt(sd_old*sd_old + (double)(data_new - data_old)*(data_new-m_newM + data_old - m_oldM)/(window_size-1));
				m_oldM = m_newM;
				sd_old = sd_new;
		}
/*
		if (j>(window_size+1+start_time) && j < (start_time+(num_samples-window_size))){
			m_newM = m_oldM - (double)data[(-window_size+j-2)*channels_per_spec+start_channel-start_channel_offset]/window_size + (double)data[(j-1)*channels_per_spec+start_channel-start_channel_offset]/window_size;
				sd_new = sqrt(sd_old*sd_old + (double)(data[(j-1)*channels_per_spec+start_channel-start_channel_offset] - data[(-window_size+j-2)*channels_per_spec+start_channel-start_channel_offset])*(data[(j-1)*channels_per_spec+start_channel-start_channel_offset]-m_newM + data[(j-window_size-2)*channels_per_spec+start_channel-start_channel_offset] - m_oldM)/(window_size-1));
				m_oldM = m_newM;
				sd_old = sd_new;
		}
*/
		//printf("sd %.6f \n", sd_new);
		//printf("mean %.6f \n", m_newM);
		//printf("number %d\n\n", data[j*channels_per_spec+start_channel-start_channel_offset]);
		avg[(j-start_time)*channels_per_spec+start_channel] = m_newM;
		sd[(j-start_time)*channels_per_spec+start_channel] = sd_new;
		threshold_top = m_newM+sd_new*3;
		threshold_bottom = m_newM-sd_new*3;
		if (data[j*channels_per_spec+start_channel-start_channel_offset] > threshold_top || data[j*channels_per_spec+start_channel-start_channel_offset] < threshold_bottom ){//&& (j > window_size+1+start_time)){
			counter++;
			#if REPLACE_NARROW_ZERO==TRUE
				for(int n = start_channel; n < stop_channel; n++){
					data[j*channels_per_spec+n] = 0;
				}
			#else
				for(int n = start_channel; n < stop_channel; n++){
	//				printf("3\n");
					random = gasdev(&seed);
					while(random>1 || random<-1){
						random = gasdev(&seed);
					}
					//printf("%d  %f\n", j, sd[(j-start_time)*channels_per_spec+n]);
					data[j*channels_per_spec+n] = (short int)(avg[(j-start_time)*channels_per_spec+start_channel]+sd[(j-start_time)*channels_per_spec+start_channel]*random);
				}
			#endif	
			}
		j++;
	}
	return counter;
}

static int myCompare (const void * a, const void * b)
{
    return strcmp (*(const char **) a, *(const char **) b);
}
 
void sort(const char *arr[], int n)
{
    qsort (arr, n, sizeof (const char *), myCompare);
}



			








