#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <dirent.h>
#include "crab_databuf.h"
#include "filterbank.h"
#include "rand_num_generator.h"

#define TRUE 1
#define FALSE 0 
#define RAW_DATA FALSE				//as opposed to filterbank format data
#define RAW_TO_FIL FALSE              		//if TRUE, turns raw data into filterbank format by adding a header
#define SINGLE_FILE TRUE
#define ZERO_NON_BANDPASS TRUE
#define REMOVE_BAD_CHANNEL TRUE
#define REPLACE_BAD_ZERO FALSE
#define REMOVE_NARROWBAND_RFI TRUE
#define REPLACE_NARROW_ZERO FALSE
#define REMOVE_BROADBAND_RFI TRUE		//if TRUE, REMOVE_NARROWBAND_RFI must be TRUE as well
#define REPLACE_BROAD_ZERO FALSE
#define HIGHPASS_FILTER FALSE

#define START_CHANNEL 1800			//1800 if spectrum contains 4096 channels, 300 (240) if 1536 channels
#define STOP_CHANNEL 2730			//2730 if spectrum contains 4096 channels, 1230 (1264) if 1536 channes
#define WINDOW_SIZE 1001
#define CHANNELS_PER_SPEC 4096			//1536

void remove_bad_channel(short int data[], int data_len, int start_channel, int stop_channel);

int calc_stats_broad(double data[], short int replace[], double avg[], double sd[], int num_samples, int start_channel_offset, int start_channel, int stop_channel, int window_size, int channels_per_spec, long start_time);

int calc_stats_narrow(short int data[], short int replace[], double avg[], double sd[], int data_len, int start_channel_offset, int start_channel, int stop_channel, int window_size, int channels_per_spec, long start_time);

char* create_name(char *src_1, const char *src_2);

int main(int argc, char *argv[]){
	char temp[255], fileName[255], out_file[255], header[363], str[10], c;
	char clean[] = "clean_";
	char fil[] = ".fil";
	int num_files=0, counter=0, check, header_len;
	long long f_size, data_size, sum_channels = 0;
	DIR *p;
	struct dirent *pp;
	FILE * fp;
	FILE * out;

	p = opendir (argv[1]);
	if (p != NULL)
	{
		while ((pp = readdir (p))!=NULL) {
		int length = strlen(pp->d_name);
			#if RAW_DATA == TRUE
				if (strncmp(pp->d_name, "data", 4) == 0 && strncmp(pp->d_name+length-4, ".fil", 4) != 0) {
					num_files++;
				}
			#endif
			#if RAW_DATA == FALSE
				if (strncmp(pp->d_name+length-4, ".fil", 4) == 0) {
					num_files++;
				}
			#endif
		}
	}
	rewinddir(p);
	printf("Number of files: %d\n", num_files);
	char files[num_files][255];

	if (p != NULL)
	{
		while ((pp = readdir (p))!=NULL) {
		int length = strlen(pp->d_name);   
			#if RAW_DATA == TRUE
	      		if (strncmp(pp->d_name, "data", 4) == 0 && strncmp(pp->d_name+length-4, ".fil", 4) != 0) {
				strncpy(fileName, pp->d_name, 254);
        			fileName[254] = '\0';
				strcpy(files[counter], fileName);
				counter++;
	      		}
			#endif
			#if RAW_DATA == FALSE
	      		if (strncmp(pp->d_name+length-4, ".fil", 4) == 0) {
				strncpy(fileName, pp->d_name, 254);
        			fileName[254] = '\0';
				strcpy(files[counter], fileName);
				counter++;
	      		}
			#endif
	   	}
	    	(void) closedir (p);
	}

	for(int i=0;i< num_files;i++)
	for(int j=i+1;j< num_files;j++){
		 if(strcmp(files[i],files[j])>0){
		    strcpy(temp,files[i]);
		    strcpy(files[i],files[j]);
		    strcpy(files[j],temp);
		 }
	}
	
	#if SINGLE_FILE == TRUE
		strcpy(out_file, create_name(clean,"data"));
		#if RAW_TO_FIL == TRUE
			strcpy(out_file, create_name(out_file, fil));
			WriteHeader(out_file);
		#endif
		out = fopen(out_file, "a+");
		if (out == NULL){
			printf("Could not open out file");
		}
	#endif
	for(int i = 0; i<num_files;i++){
		check = FALSE;
		fp = fopen(files[i], "rb");
		printf("Input filename: %s\n", files[i]);
		fseek(fp, 0L, SEEK_END);
		f_size = ftell(fp);
		rewind(fp);
		if (f_size != 0){//(f_size <= (246345728*8) && f_size != 0){  //num data points per minute of data
		printf("File size: %lld\n", f_size);
			#if SINGLE_FILE == FALSE
				strcpy(out_file, create_name(clean, files[i]));
				#if RAW_TO_FIL == TRUE
					strcpy(out_file, create_name(out_file,fil));
					WriteHeader(out_file);	
				#endif
				out = fopen(out_file, "a+");
			#endif
			fseek(out, 0L, SEEK_END);
			header_len = ftell(out);
			#if RAW_DATA == FALSE
				while (fgets (&c, 1, fp)!=NULL) {
					fread (&c, 1, 1, fp);
					if(c == 'H'){
						fseek(fp, -1, SEEK_CUR);
						fread(str, 10, 1, fp);
						if(strstr(str, "HEADER_END") != NULL){
							break;
						}
					}
				   }
				header_len = ftell(fp);
				rewind(fp);
				if (ftell(out) < 100 && f_size <= (246345728*8)){
					fread(header, sizeof(char), header_len, fp);
					fwrite(header, sizeof(char), header_len, out);
				}
			#endif
			printf("Header size = %d bytes\n", header_len);
			fseek(fp, header_len, SEEK_SET);
			data_size = f_size-header_len;
			printf("Data size =  %lld bytes\n", data_size);
			if (data_size != 0){
				check = TRUE;
				printf("Writing data to: %s\n", out_file);
				long n_data_points = data_size/sizeof(short int);
				long num_time_samples = n_data_points/CHANNELS_PER_SPEC;
				printf("Number of time samples: %ld\n", num_time_samples);

				short int *data = (short int*)malloc(sizeof(short int)*n_data_points);
				#if REMOVE_BROADBAND_RFI == TRUE
					double *time_sum = (double*)malloc(sizeof(double)*num_time_samples);
					if (!data ||!time_sum){
						printf("Malloc Failed");
					}
				#endif
				#if REMOVE_NARROWBAND_RFI == TRUE
					double *channel_avg = (double*)malloc(sizeof(double)*n_data_points);
					double *channel_SD = (double*)malloc(sizeof(double)*n_data_points);
					if (!data || !channel_avg || !channel_SD){
						printf("Malloc Failed");
					}
				#endif

				if (!data){
					printf("Malloc Failed");
				}

				int size = fread(data, sizeof(short int), n_data_points, fp);

				#if ZERO_NON_BANDPASS == TRUE	
				for(int j = 0; j < num_time_samples; j++){
					for(int k = 0; k< CHANNELS_PER_SPEC; k++){
						if (k < START_CHANNEL || k >= STOP_CHANNEL){
							data[CHANNELS_PER_SPEC*j+k] = 0;
						}
					}
				}
				#endif

				#if REMOVE_BAD_CHANNEL == TRUE
					remove_bad_channel(data, n_data_points, START_CHANNEL, STOP_CHANNEL);
				#endif
				
				#if REMOVE_NARROWBAND_RFI == TRUE
					counter = 0;
					printf("Removing narrowband RFI\n");
					for(int n = START_CHANNEL; n < STOP_CHANNEL; n++){
						counter += calc_stats_narrow(data, data, channel_avg, channel_SD, num_time_samples, 0, n, n+1, WINDOW_SIZE, CHANNELS_PER_SPEC,0);
					}
					for(int n = START_CHANNEL; n < STOP_CHANNEL; n++){
						counter += calc_stats_narrow(data, data, channel_avg, channel_SD, num_time_samples, 0, n, n+1, WINDOW_SIZE, CHANNELS_PER_SPEC,0);
					}
					printf("Number of narrowband samples removed: %d\n", counter);
				#endif
				
				#if REMOVE_BROADBAND_RFI == TRUE
					counter = 0;
					printf("Removing broadband RFI\n");
					for(int j = 0; j < num_time_samples; j++){
						for(int k = START_CHANNEL; k< STOP_CHANNEL; k++){
							//printf("DATA VALUE: %d\n", I[num_channel_per_spec*j+k]);
							//printf("AVG VALUE %d\n\n", (short int)channel_avg[num_channel_per_spec*j+k]);
							sum_channels += data[CHANNELS_PER_SPEC*j+k];
						}
						time_sum[j] = sum_channels;
			//			printf("TIME_SUM: %f\n\n", time_sum[j]);
						sum_channels = 0;
					}
					counter = calc_stats_broad(time_sum, data, channel_avg, channel_SD, num_time_samples, START_CHANNEL, START_CHANNEL, STOP_CHANNEL, WINDOW_SIZE, 1, 0);
					printf("Number of broadband samples removed: %d\n",counter);
				#endif

				#if HIGHPASS_FILTER == TRUE
					for(int j = 0; j < num_samples; j++){
						for(int k = 0; k< CHANNELS_PER_SPEC; k++){
						//printf("DATA VALUE: %d\n", I[num_channel_per_spec*j+k]);
						//printf("AVG VALUE %d\n\n", (short int)channel_avg[num_channel_per_spec*j+k]);
						data[CHANNELS_PER_SPEC*j+k]-= (short int)channel_avg[CHANNELS_PER_SPEC*j+k];
						}
					}
				#endif

				fwrite(data,sizeof(short int),n_data_points,out);
				free(data);
				#if REMOVE_NARROWBAND_RFI == TRUE
					free(channel_avg);
					free(channel_SD);
				#endif
				#if REMOVE_BROADBAND_RFI == TRUE
					free(time_sum);
				#endif
				//fread(data, 1, data_size, fp);
				//fwrite(data, 1, data_size, out);
				//free(data);
			}
			#if SINGLE_FILE == FALSE
				if(check == FALSE){
					remove(out_file);
				}
				fclose(out);
			#endif
		}
		if(check == FALSE){
			printf("Did not write data to file\n");
		}
		fclose(fp);
		printf("\n");
	}
	#if SINGLE_FILE == TRUE
		fclose(out);
	#endif
	return 0;
}

char* create_name(char *src_1, const char *src_2){
	size_t src_2_len = strlen(src_2);
	size_t src_1_len = strlen(src_1);
	size_t i;
	char *output = (char *)malloc(src_2_len+src_1_len);
	for (i = 0 ; i < src_1_len && src_1[i] != '\0' ; i++)
		output[i] = src_1[i];
	for(i = src_1_len ; i < (src_1_len+src_2_len) && src_2[(i-src_1_len)] != '\0'; i++)
		output[i] = src_2[(i-src_1_len)];
	output[i] = '\0';

	return output;
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
			m_newM = m_oldM - (double)data[(-window_size+j-1)*channels_per_spec+start_channel-start_channel_offset]/window_size + (double)data[(j)*channels_per_spec+start_channel-start_channel_offset]/window_size;
				sd_new = sqrt(sd_old*sd_old + (double)(data[(j)*channels_per_spec+start_channel-start_channel_offset] - data[(-window_size+j-1)*channels_per_spec+start_channel-start_channel_offset])*(data[(j)*channels_per_spec+start_channel-start_channel_offset]-m_newM + data[(j-window_size-1)*channels_per_spec+start_channel-start_channel_offset] - m_oldM)/(window_size-1));
				if(isnan(sd_new)){
					sd_new = sd_old;
				}
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
		if (data[j*channels_per_spec+start_channel-start_channel_offset] > threshold_top || data[j*channels_per_spec+start_channel-start_channel_offset] < threshold_bottom ){//&& j>(window_size+start_time)){
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
		if (j>(window_size+1+start_time) && j < (start_time+(num_samples-window_size))){
			m_newM = m_oldM - (double)data[(-window_size+j-2)*channels_per_spec+start_channel-start_channel_offset]/window_size + (double)data[(j-1)*channels_per_spec+start_channel-start_channel_offset]/window_size;
				sd_new = sqrt(sd_old*sd_old + (double)(data[(j-1)*channels_per_spec+start_channel-start_channel_offset] - data[(-window_size+j-2)*channels_per_spec+start_channel-start_channel_offset])*(data[(j-1)*channels_per_spec+start_channel-start_channel_offset]-m_newM + data[(j-window_size-2)*channels_per_spec+start_channel-start_channel_offset] - m_oldM)/(window_size-1));
				if(isnan(sd_new)){
					sd_new = sd_old;
				}
				m_oldM = m_newM;
				sd_old = sd_new;
		}
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
