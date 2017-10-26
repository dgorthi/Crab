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
#include "acq_udp.h"

#define TRUE 1
#define FALSE 0 
#define RAW_DATA TRUE				//as opposed to filterbank format data
#define RAW_TO_FIL TRUE             		//if TRUE, turns raw data into filterbank format by adding a header
#define SINGLE_FILE TRUE			//if TRUE, combines all input files to into a single output file

#define ZERO_NON_BANDPASS TRUE
#define REMOVE_BAD_CHANNEL TRUE
#define REPLACE_BAD_ZERO FALSE

#define REMOVE_NARROWBAND_RFI TRUE
#define REPLACE_NARROW_ZERO FALSE
#define PRINT_NARROW_STATS FALSE		//avg, sd, threshold, replacements, etc
#define PRINT_NARROW_STATS1 FALSE		//which time samples have bad values

#define REMOVE_BROADBAND_RFI TRUE		//if TRUE, REMOVE_NARROWBAND_RFI must be TRUE as well
#define REPLACE_BROAD_ZERO FALSE
#define PRINT_BROAD_STATS FALSE			//avg, sd, threshold, replacements, etc.
#define PRINT_BROAD_STATS1 FALSE		//which time samples have bad values

#define BOTH_POLS TRUE
#define XX_POL FALSE
#define YY_POL FALSE

#define START_CHANNEL 240			//1800 if spectrum contains 4096 channels, 300 (240) if 1536 channels
#define STOP_CHANNEL 1264			//2730 if spectrum contains 4096 channels, 1230 (1264) if 1536 channes 2270-2730
#define WINDOW_SIZE 1000
#define CHANNELS_PER_SPEC 1536			//1536 or 4096


char* create_name(char *src_1, const char *src_2);

void remove_bad_channel(short int data[], int data_len, int start_channel, int stop_channel);

int remove_broadband_RFI(double input[], short int output[], double avg[], double sd[], long start_time, long num_samples, int channel, int channels_per_samp, long seed, double thresh_limit);

int remove_narrowband_RFI(short int input[], short int output[], double avg[], double sd[], long start_time, long num_samples, int channel, int channels_per_samp, long seed, double thresh_limit);

int main(int argc, char *argv[]){
	/*This first part of the main function handles the file input of the program. The user must input the directory the data files 
	they wish to convert to filterbank and/or remove RFI are in  when they call the program (remove_RFI.c /directory_of_files). If the
	files are al=ready in filterbank then the user should set the RAW_DATA macro to false. Otherwise, if the data is in binary format, 
	the user should set this macro to true. If the user would additionally like to convert the raw data to filterbank format, they should
	set the RAW_TO_FIL parameter to true. Lastly, if the user would like to combine multiple input files into a single output file, they
	should set the SINGLE_FILE macro to be true. 
	NOTE: raw data files must have the first four characters of the file name be data, and filterbank files must have the last four characters
	be .fil
	*/

        char temp[255], fileName[255], out_file[255], header[363], str[10], c;
	char *clean = (char *)malloc(sizeof(char)*255);
	#if BOTH_POLS == TRUE
        	clean = create_name(argv[1],"reduced/clean1_");
	#endif
	#if XX_POL == TRUE
        	clean = create_name(argv[1],"reduced/xx_clean_");
	#endif
	#if YY_POL == TRUE
        	clean = create_name(argv[1],"reduced/yy_clean_");
	#endif
	//fprintf(stderr,"%s\n",clean);
	char fil[] = ".fil";
	int sample_size = UDP_DATA*8, num_fft_channels = 4096, sampling_rate = 2100000000, acc_len = pow(2,11), num_files=0, check, header_len, no_cpy_buf, size, counter = 0;
	short int  xx, yy, xy_re, xy_im;
	double random, sampling_time = 2*num_fft_channels*acc_len/(double)sampling_rate;
	long n_data_points, num_time_samples, seed = 5;
	long long f_size, data_size, sum_channels = 0;
	DIR *p;
	struct dirent *pp;
	FILE * fp;
	FILE * out;
	/*Opens file directory and counts number of data files in directory*/
	p = opendir(argv[1]);
	fprintf(stderr,"Opened the directory!\n");
	if (p != NULL)
	{
		while ((pp = readdir (p))!=NULL) {
		int length = strlen(pp->d_name);
		fprintf(stderr,"Opening: %s..\n",pp->d_name);
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
	printf("Rewinding!\n");
	rewinddir(p);
	fprintf(stderr,"Number of data files in directory: %d\n\n", num_files);
	char files[num_files][255];

	/*Creates an array of the filenames for future access*/
	if (p != NULL)
	{
		while ((pp = readdir (p))!=NULL) {
		int length = strlen(pp->d_name);   
			#if RAW_DATA == TRUE
	      		if (strncmp(pp->d_name, "data", 4) == 0 && strncmp(pp->d_name+length-4, ".fil", 4) != 0 && strncmp(pp->d_name+length-4, ".txt", 4) != 0) {
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

	/*Reorders the files in the array by time (DIR randomly accesses the files in the directory so the original array of 
	filenames is all out of order*/
	for(int i=0;i< num_files;i++)
	for(int j=i+1;j< num_files;j++){
		 if(strcmp(files[i],files[j])>0){
		    strcpy(temp,files[i]);
		    strcpy(files[i],files[j]);
		    strcpy(files[j],temp);
		 }
	}
	
	/*If user wants to combine multiple input files into a single file the file is created now. If the user is converting
	raw data to filterbank format, a header is added to the file. Otherwise the header of the first filterbnak file is
	used as the header for the combined file.
	NOTE: If converting from raw data, the time information in the header will correspond to the time this program is run
	It needs to be fixed manually to match the time the data was actually collected.*/
	#if SINGLE_FILE == TRUE
		strcpy(out_file, create_name(clean,"data"));
		//strcpy(out_file, create_name(argv[1], out_file));
		#if RAW_TO_FIL == TRUE
			strcpy(out_file, create_name(out_file, fil));
			fprintf(stderr,"OUTPUT FILE NAME: %s\n\n",out_file);
			WriteHeader(out_file);
		#endif
		out = fopen(out_file, "a+");
		if (out == NULL){
			printf("Failed to open output file\n\n\nPOSSIBLY DUE TO DIRECTORY PERMISSIONS!!!!\nCHMOD 777!!!!!\n");
		}
	#endif

	/*This part of the code removes RFI from the input data files and creates the filterbank file(s). The RFI removal
	uses a basic threshold test for both broadband and narrowband RFI. 
	NOTE: Uses Welford's method for calculation of standard deviation.*/ 
	for(int i = 0; i<num_files;i++){
		check = FALSE;
		//Opens ith file in directory
		fp = fopen(create_name(argv[1],files[i]), "rb");
		if(fp == NULL){
			printf("Failed to open input file\n\n\nPOSSIBLY DUE TO DIRECTORY PERMISSIONS!!!!\nCHMOD 777!!!!!\n");
		}
		printf("INPUT FILENAME: %s\n", create_name(argv[1],files[i]));
		fseek(fp, 0L, SEEK_END);
		f_size = ftell(fp);
		rewind(fp);
		//Checks to see if empty
		if (f_size != 0){//(f_size <= (246345728*8) && f_size != 0){  //num data points per minute of data
			printf("File size: %lld\n", f_size);
			#if SINGLE_FILE == FALSE
				strcpy(out_file, create_name(clean, files[i]));
				#if RAW_TO_FIL == TRUE
					strcpy(out_file, create_name(out_file,fil));
					WriteHeader(out_file);	
				#endif
				out = fopen(out_file, "a+");
				printf("OUTPUT FILENAME: %s\n", out_file);
				if (out == NULL){
					printf("Failed to open output file\n\n\nPOSSIBLY DUE TO DIRECTORY PERMISSIONS!!!!\nCHMOD 777!!!!!\n");
				}
			#endif
			fseek(out, 0L, SEEK_END);
			//header_len = ftell(out);
			header_len = ftell(fp);
			//Copies original header of old filterbank file to copy to new output file
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
				no_cpy_buf = data_size/(ACC_BUFSIZE+NACC);
				//n_data_points = (data_size-NACC*no_cpy_buf)/sizeof(short int);
				n_data_points = NACC*UDP_DATA*no_cpy_buf/sizeof(short int);
				//num_time_samples = 3/8.*n_data_points/CHANNELS_PER_SPEC;
				num_time_samples = n_data_points*2*8/3/(sample_size);
				printf("Number of copy buffers: %d\n", no_cpy_buf);
				printf("Number of data points: %ld\n", n_data_points);
				printf("Number of time samples: %ld\n", num_time_samples);

				//Allocates memory for data array
				short int *data_full = (short int*)malloc(sizeof(short int)*n_data_points);
				short int *data = (short int*)malloc(sizeof(short int)*n_data_points/4);
				char *flag = (char*)malloc(sizeof(char) * NACC);  
				#if REMOVE_BROADBAND_RFI == TRUE
					double *time_sum = (double*)malloc(sizeof(double)*num_time_samples);
					if (!data ||!time_sum){
						printf("Malloc Failed, input data size possibly too large");
					}
				#endif
				#if REMOVE_NARROWBAND_RFI == TRUE
					double *channel_avg = (double*)malloc(sizeof(double)*n_data_points/4);
					double *channel_SD = (double*)malloc(sizeof(double)*n_data_points/4);
					if (!data || !channel_avg || !channel_SD){
						printf("Malloc Failed, input data size possibly too large");
					}
				#endif

				if (!data || !flag || !data_full){
					printf("Malloc Failed, input data size possibly too large");
				}
				
				//Reads data into an array
//				int size = fread(data, sizeof(short int), n_data_points, fp);
				for (int ii = 0; ii < no_cpy_buf; ii++){
//					printf("Copy buffer %d\n", ii);
/*					
					for (int jj = 0; jj < ACC_BUFSIZE/2; jj++){
						size = fread(&xx, sizeof(short int), 1, fp);
						size = fread(&yy, sizeof(short int), 1, fp);
						fseek(fp, +2, SEEK_CUR);
						//size = fread(&xy_re, sizeof(short int), 1, fp);
						//size = fread(&xy_im, sizeof(short int), 1, fp);
						data_full[ii*ACC_BUFSIZE/8+jj] = (short int)sqrt(double(xx)*double(xx)+(double)yy*(double)yy);
	//					printf("Data point: %d, Value: %d\n", jj, data[ii*ACC_BUFSIZE/8+jj]);
	//					printf("XX: %d, YY: %d, XY_RE: %d, XY_IM: %d, \n\n", xx, yy,xy_re, xy_im);
					}
*/					
					size = fread((data_full+ii*ACC_BUFSIZE/2), sizeof(short int), ACC_BUFSIZE/2, fp);
					size = fread(flag, sizeof(char), NACC, fp);
					for (int jj = 0; jj < NACC; jj++){
						if (flag[jj] == 1 && ii == 0){
							printf("Packet number: %d missed, replacing with zeros\n", jj);
							for (int kk = 0; kk < 512*4; kk++){
                                                                data_full[ii*ACC_BUFSIZE/2+512*4*jj+kk] = 0;
                                                        }
						}
						if (flag[jj] == 1 && ii != 0){
							printf("Packet number: %d missed, replacing with guassian noise\n", jj);
							for (int kk = 0; kk < 512*4; kk++){ 
								/*while(random>.2 || random<-.2){
									random = gasdev(&seed);
								}*/
	//							data_full[ii*ACC_BUFSIZE+jj*UDP_DATA+kk] = data_full[(ii-1)*ACC_BUFSIZE+jj*UDP_DATA+kk]*random;
								data_full[ii*ACC_BUFSIZE/2+512*4*jj+kk] = data_full[(ii-1)*ACC_BUFSIZE/2+jj*512*4+kk]*random;
							}
						}
					}
				}

				for (int i = 0; i < n_data_points/4-1; i++){
					#if BOTH_POLS == TRUE
						data[i] = sqrt(data_full[4*i]*data_full[4*i]+data_full[4*i+1]*data_full[4*i+1]);
					#endif
					#if XX_POL == TRUE
						data[i] = data_full[4*i];
					#endif
					#if YY_POL == TRUE
						data[i] = data_full[4*i+1];
					#endif
				//	printf("%d\n", data[i]);
				}
/*
				for (int i = 10000; i < 11000; i=i+4){
					fprintf(stderr,"Data point %d-%d: XX %d, YY %d, XY_RE %d, XY_IM %d\n", i,i+3, data_full[i], data_full[i+1], data_full[i+2], data_full[i+3]);
					fprintf(stderr,"Power %d\n\n", data[i/4]);
				}
*/

				printf("Setting values outside bandpass to zero...\n");
				#if ZERO_NON_BANDPASS == TRUE	
				for(int j = 0; j < num_time_samples; j++){
					for(int k = 0; k< CHANNELS_PER_SPEC; k++){
						if (k < START_CHANNEL || k >= STOP_CHANNEL){
							data[CHANNELS_PER_SPEC*j+k] = 0;
						}
					}
				}
				#endif
				printf("Setting values outside bandpass to zero done\n");


				printf("Removing bad channels...\n");
				#if REMOVE_BAD_CHANNEL == TRUE
					remove_bad_channel(data, n_data_points/4, START_CHANNEL, STOP_CHANNEL);
				#endif
				printf("Removing bad channels done\n");

				free(data_full);
				free(flag);

			/*
			RFI removal is done by a threshold test. For each channel, if the current time sample of interest
			is > WINDOW_SIZE+start_time sample  number of time samples, the average and standard deviation for the threshold test is 
			calculated from the preceeding WINDOW_SIZE number of time samples. If the current time sample of interest
			is within the first WINDOW_SIZE+start_time number of time samples, then the statistics for that channel at
			that time sample for the threshold test are calculated from the first WINDOW_SIZE number of time samples. 
			If the current time sample of interest is within the last WINDOW_SIZE number of time samples, then the
			statistics for that channel at that time sample for the threshold test are calculated from the last
			WINDOW_SIZE number of time samples.
			When the statistics have been calculated the threshold limit is set at +/-3 standard deviations
			from the average, and the power of the channel at the current time sample is compared to the threshold 
			limit. If it lies outside of the threshold limit, then the power is replaced either by 0 or by 
			gaussian noise about the average value used in the calculation of the threshold limit depending on the
			macro setting. 

			The arguments for the function are (in order) the following:
			
			int remove_RFI(double input[], short int output[], double avg[], double sd[], long start_time, long num_samples, int channel, int channels_per_samp);

			double input[]: The array that the data is read from (RFI not yet removed)
			int output[]: The array that the data is output to (RFI removed)
			double avg[]: The array which holds the average value of each channel at each time sample as calculated according 
			to the method above
			double sd[]: The array which holds the standard deviation of each channel at each time sample as calculated according
			to the method above
			long start_time: Start time sample to remove RFI from
			long num_samples: Number of time samples after start_time to remove RFI from
			int channel: Channel to remove RFI from 
			int channels_per_samp: The number of channels each time sample of the data array contains
			double thresh_limit: number of standard deviations from mean with which to set threshold limits

			NOTE: The remove_RFI function can only handle one channel at a time
			*/  
				#if REMOVE_NARROWBAND_RFI == TRUE
					/*
					For narrowband RFI, when a value of the input array exceeds the threshold limit for that value
					the value stored in its memory location of that value is overridden  with gaussian noise or 0. Thus the 
					input data array is the same as the output data array- the input having RFI, the output with RFI
					removed.  
					The channels_per_samp is thus CHANNELS_PER_SPEC, as each time sample of the input array has CHANNELS_PER_SPEC
					channels (contrast to broadband RFI removal).
					The for loop iterates through every channel in our bandpass and removes the RFI for that channel. A counter
					of the number of time samples replaced of that channel is then returned, and the sum of these counters (total
					number of data points that had narrowband RFI) is printed.
					*/ 
					counter = 0;
					printf("Removing narrowband RFI\n");
					for(int n = START_CHANNEL; n < STOP_CHANNEL; n++){
						counter += remove_narrowband_RFI(data, data, channel_avg, channel_SD, 0, num_time_samples, n, CHANNELS_PER_SPEC, seed, 3);
	//					counter += calc_stats_narrow(data, data, channel_avg, channel_SD, num_time_samples, 0, n, n+1, WINDOW_SIZE, CHANNELS_PER_SPEC,0);
						//printf("Channel %d\n",n);
					}
					printf("Number of narrowband samples removed: %d\n", counter);
				#endif


				#if REMOVE_BROADBAND_RFI == TRUE
					/*
					Broadband RFI removal works much like narrowband RFI removal except that our "channel" is really the sum of
					all our channels across the band for a given time sample. Thus our input data array is an array with the total
					power of each time sample at each index, and our CHANNELS_PER_SPEC is 1.
					Our output array, however, is the same as the output array for the narrowband RFI removal. If any "channel" exceeds 
					the threshold limit for that "channel" at the given time sample, then every channel of that time sample is replaced
					with either 0 or gaussian noise.
					
					NOTE: If REMOVE_BROADBAND_RFI is set to TRUE, then so must REMOVE_NARROWBAND_RFI, as the statistics used in 
					calculating the gaussian noise to replace the bad values with our calculated in the REMOVE_NARROWBAND_RFI
					part.  
					*/ 
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
					counter = remove_broadband_RFI(time_sum, data, channel_avg, channel_SD, 0, num_time_samples, 0, 1, seed, 3);
//					counter = calc_stats_broad(time_sum, data, channel_avg, channel_SD, num_time_samples, START_CHANNEL, START_CHANNEL, STOP_CHANNEL, WINDOW_SIZE, 1, 0);
					printf("Number of broadband samples removed: %d\n",counter);
				#endif
				
				printf("Writing cleaned data to file %s...\n", out_file);
				fwrite(data,sizeof(short int),n_data_points/4,out);
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
			printf("Did not write data to file. Data file possibly empty.\n");
		}
		fclose(fp);
		printf("\n");
	}
	#if SINGLE_FILE == TRUE
		fclose(out);
	#endif
	return 0;
}

/*The original strcpy function replaces the string in the first argument with the combined string, whereas
this function produces a new string that is a combination of both arugments, not changing either of the
arguments*/
char* create_name(char *src_1, const char *src_2){
	size_t src_2_len = strlen(src_2);
	size_t src_1_len = strlen(src_1);
	size_t i;
	char *output = (char *)malloc(src_2_len+src_1_len+2);
	for (i = 0 ; i < src_1_len && src_1[i] != '\0' ; i++)
		output[i] = src_1[i];
	for(i = src_1_len ; i < (src_1_len+src_2_len) && src_2[(i-src_1_len)] != '\0'; i++)
		output[i] = src_2[(i-src_1_len)];
	output[i] = '\0';
	return output;
}

/*Removes bad channels by comparing value of each channel within a single time sample to a threshold
set by the avg and the sd of the channels at that time sample*/
void remove_bad_channel(short int data[], int data_len, int start_channel, int stop_channel){
	int num_channel = 1536;
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
		thresh_top = m_newM + 2*sd_new;
		thresh_bottom = m_newM-2*sd_new;
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
					/*while(random>1 || random<-1){
						random = gasdev(&seed);
					}*/
					data[j*num_channel+n] = (short int)(m_newM+sd_new*random);
				}	
			}
		#endif	
	}
}

/*Removes narrowband RFI by simple threshold test*/
int remove_narrowband_RFI(short int input[], short int output[], double avg[], double sd[], long start_time, long num_samples, int channel, int channels_per_samp, long seed, double thresh_limit){
	/*
	If time sample is < start_time+WINDOW_SIZE, the average and standard deviation are calculated according to the algorithm described in 
	http://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
	If the time sample is >= start_time+WINDOW_SIZE, a rolling average and standard deviation according to the algorithm outlined at 
	http://jonisalonen.com/2014/efficient-and-accurate-rolling-standard-deviation/ 
	*/
	int counter = 0;
	long next_index, current_index, remove_index;
	double old_sample, current_sample, next_sample, oldM, M, S, V, SD, thresh_top, thresh_bottom, random;

	/*Set aside memory for values used in calculation of statistics*/ 
	double *window = (double*)malloc(sizeof(double)*WINDOW_SIZE);
	if (!window){
		printf("Malloc Failed, possibly out of memory");
	}


	for(long j = start_time; j < (num_samples-start_time-1); j++){
		if (j == start_time){
			/*Calculate stats for first WINDOW_SIZE time samples*/
			M = 0;
			S = 0;
			for(int k = 1; k <= WINDOW_SIZE; k++){
				current_index = (k-1)*channels_per_samp+channel;
				current_sample = (double)input[current_index];
				oldM = M;
				M = M + (current_sample-M)/k;
				S = S + (current_sample-M)*(current_sample-oldM);
			}
			V = S/(WINDOW_SIZE-1);
			SD = sqrt(V);
			thresh_top = M + thresh_limit*S;
			thresh_bottom = M - thresh_limit*S;

			/*Compare first thousand time samples to threshold limits and replace with either 0 or noise if outside*/
			for(int k = 0; k < WINDOW_SIZE; k++){
				current_index = k*channels_per_samp+channel;
				current_sample = (double)input[current_index];
				window[k] = current_sample;
				if(current_sample > thresh_top || current_sample < thresh_bottom){
					#if PRINT_NARROW_STATS1 == TRUE
						printf("******************************************************************\n");
						printf("HITTTT!!!!!!!!\n");
						printf("Current time sample:%ld, and channel: %d\n", j, channel);
						printf("*******************************************************************\n");
					#endif
					#if PRINT_NARROW_STATS == TRUE
						printf("HIT!!!!!!!!!!\n");
						printf("Current time sample:%d, and channel: %d\n", k, channel);
						printf("Original mean: %f\n", M);
						printf("Original standard deviation: %f\n", SD);
						printf("Input value: %d\n", input[current_index]);
						printf("Original window value: %f", window[k]);
					#endif
					/*Increase counter for RFI detection*/
					counter++;
					
					/*Replace bad stat window value with noise*/
					random = gasdev(&seed);
					/*while(random>1 || random<(-1)){
						random = gasdev(&seed);
					}*/
					window[k] = M + SD*random;
					#if REPLACE_NARROW_ZERO == TRUE
							output[current_index] = 0;
					#endif
					#if REPLACE_NARROW_ZERO == FALSE
							output[current_index] = (short int)window[k];
					#endif
					#if PRINT_NARROW_STATS == TRUE
						printf("Output value: %d\n", output[current_index]);
						printf("Final window value: %f\n\n ", window[k]);
					#endif
				}
			}

			/*Recalculate stats for cleaned first thousand time samples (found in window array)*/
			for(int k = 1; k <= WINDOW_SIZE; k++){
				current_sample = (double)window[k-1];
				oldM = M;
				M = M + (current_sample-M)/k;
				S = S + (current_sample-M)*(current_sample-oldM);
			}
			V = S/(WINDOW_SIZE-1);
			SD = sqrt(V);

			/*Update statistics array*/
			for(int k = 0; k < WINDOW_SIZE; k++){
				current_index = k*channels_per_samp+channel;
				avg[current_index] = M;
				sd[current_index] = SD;
			}
			oldM = M;
			
			#if PRINT_NARROW_STATS == TRUE
				for(int k = 0; k < WINDOW_SIZE; k++){
					current_index = k*channels_per_samp+channel;
					printf("Current time sample:%d, and channel: %d\n", k, channel);
					printf("Final mean value: %f\n", avg[current_index]);
					printf("Final standard deviation value: %f\n\n", sd[current_index]);
				}
			#endif
		}

		/*Begin calculation of rolling avg and sd*/
		if (j >= start_time + WINDOW_SIZE && j < num_samples - WINDOW_SIZE){
			#if PRINT_NARROW_STATS == TRUE
				printf("Current time sample:%ld, and channel: %d\n", j, channel);
			#endif

			/*Determine old, current, and next samples*/
			current_index = j*channels_per_samp+channel;
			next_index = (j+1)*channels_per_samp+channel;
			old_sample = window[0];
			current_sample = window[WINDOW_SIZE-1];
			next_sample = input[next_index];
			
			/*Update standard deviation and mean*/
			oldM = M;
			M = oldM + (current_sample-old_sample)/WINDOW_SIZE;
			V = V + (current_sample-old_sample)*(current_sample-M+old_sample-oldM)/(WINDOW_SIZE-1);
			SD = sqrt(V);

			/*Update statistics arrays*/
			avg[current_index] = M;
			sd[current_index] = SD;

   			/*Shift all elements in window by one*/ 
			memmove(&window[0], &window[1], (WINDOW_SIZE)*sizeof(double));
			window[WINDOW_SIZE-1] = next_sample;
			

			/*Set threshold limits*/
			thresh_top = M + thresh_limit*SD;
			thresh_bottom = M - thresh_limit*SD;

			#if PRINT_NARROW_STATS == TRUE
				printf("Mean: %f\n", M);
				printf("Standard deviation: %f\n", SD);
				printf("Input value: %d\n", input[current_index]);
				printf("Replaced input value (for stat calc): %f\n", current_sample);
				printf("Next value: %f\n", next_sample);
				printf("Threshold top: %f\n", thresh_top);
				printf("Threshold bottom: %f\n\n", thresh_bottom);
			#endif
			/*Replace bad values*/
			if(next_sample > thresh_top || next_sample < thresh_bottom){
				#if PRINT_NARROW_STATS1 == TRUE
					printf("******************************************************************\n");
					printf("HITTTT!!!!!!!!\n");
					printf("Current time sample:%ld, and channel: %d\n", j, channel);
					printf("*******************************************************************\n");
				#endif
				/*Increase counter for RFI detection*/
				counter++;
				
				/*Replace bad stat window value with noise*/	
				random = gasdev(&seed);
				/*while(random>1 || random<(-1)){
					random = gasdev(&seed);
				}*/
				window[WINDOW_SIZE-1] = M + SD*random;

				#if REPLACE_NARROW_ZERO == TRUE
					output[next_index] = 0;
				#endif
				#if REPLACE_NARROW_ZERO == FALSE
					output[next_index] = (short int)window[WINDOW_SIZE-1];
				#endif

				#if PRINT_NARROW_STATS == TRUE
					printf("Output value: %d\n\n", output[current_index]);
				#endif
			}
			/*
			for(int k = 0; k < WINDOW_SIZE; k++){
				printf("%f, ", window[k]);
			}
				printf("_________________________________\n");*/
		}
		
		/*For last WINDOW_SIZE samples, use stats calculated from the last WINDOW_SIZE samples before sample num_samples-WINDOW_SIZE*/
		if (j >= WINDOW_SIZE){
			current_index = j*channels_per_samp+channel;
			next_index = (j+1)*channels_per_samp+channel;
			next_sample = input[next_index];

			avg[current_index] = M;
			sd[current_index] = SD;

			if(next_sample > thresh_top || next_sample < thresh_bottom){
				#if PRINT_NARROW_STATS1 == TRUE
					printf("******************************************************************\n");
					printf("HITTTT!!!!!!!!\n");
					printf("Current time sample:%ld, and channel: %d\n", j, channel);
					printf("*******************************************************************\n");
				#endif

				/*Increase counter for RFI detection*/
				counter++;
				
				/*Replace bad stat window value with noise*/	
				random = gasdev(&seed);
				/*while(random>1 || random<(-1)){
					random = gasdev(&seed);
				}*/
				window[WINDOW_SIZE-1] = M + SD*random;

				#if REPLACE_NARROW_ZERO == TRUE
					output[next_index] = 0;
				#endif
				#if REPLACE_NARROW_ZERO == FALSE
					output[next_index] = (short int)window[WINDOW_SIZE-1];
				#endif
			}
		}
		
	}

	free(window);
	return counter;

}
					
/*Removes broadband RFI by simple threshold test*/
int remove_broadband_RFI(double input[], short int output[], double avg[], double sd[], long start_time, long num_samples, int channel, int channels_per_samp, long seed, double thresh_limit){
	/*
	If time sample is < start_time+WINDOW_SIZE, the average and standard deviation are calculated according to the algorithm described in 
	http://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
	If the time sample is >= start_time+WINDOW_SIZE, a rolling average and standard deviation according to the algorithm outlined at 
	http://jonisalonen.com/2014/efficient-and-accurate-rolling-standard-deviation/ 
	*/
	int counter = 0;
	long next_index, current_index, remove_index;
	double old_sample, current_sample, next_sample, oldM, M, S, V, SD, thresh_top, thresh_bottom, random;

	/*Set aside memory for values used in calculation of statistics*/ 
	double *window = (double*)malloc(sizeof(double)*WINDOW_SIZE);
	if (!window){
		printf("Malloc Failed, possibly out of memory");
	}


	for(long j = start_time; j < (num_samples-start_time-1); j++){
		if (j == start_time){
			/*Calculate stats for first WINDOW_SIZE time samples*/
			M = 0;
			S = 0;
			for(int k = 1; k <= WINDOW_SIZE; k++){
				current_index = (k-1)*channels_per_samp+channel;
				current_sample = (double)input[current_index];
				oldM = M;
				M = M + (current_sample-M)/k;
				S = S + (current_sample-M)*(current_sample-oldM);
			}
			V = S/(WINDOW_SIZE-1);
			SD = sqrt(V);
			thresh_top = M + thresh_limit*S;
			thresh_bottom = M - thresh_limit*S;

			/*Compare first thousand time samples to threshold limits and replace with either 0 or noise if outside*/
			for(int k = 0; k < WINDOW_SIZE; k++){
				current_index = k*channels_per_samp+channel;
				current_sample = (double)input[current_index];
				window[k] = current_sample;
				if(current_sample > thresh_top || current_sample < thresh_bottom){
					#if PRINT_BROAD_STATS1 == TRUE
						printf("******************************************************************\n");
						printf("HITTTT!!!!!!!!\n");
						printf("Current time sample: %ld\n", j);
						printf("*******************************************************************\n");
					#endif
					#if PRINT_BROAD_STATS == TRUE
						printf("HIT!!!!!!!!!!\n");
						printf("Current time sample: %d\n", j);
						printf("Original mean: %f\n", M);
						printf("Original standard deviation: %f\n", SD);
						printf("Input value: %f\n", input[current_index]);
						printf("Original window value: %f", window[k]);
					#endif
					/*Increase counter for RFI detection*/
					counter++;
					
					/*Replace bad stat window value with noise*/
					random = gasdev(&seed);
					/*while(random>1 || random<(-1)){
						random = gasdev(&seed);
					}*/
					window[k] = M + SD*random;

					#if REPLACE_BROAD_ZERO == TRUE
						for(int l = START_CHANNEL; l < STOP_CHANNEL; l++){
							next_index = (j+1)*CHANNELS_PER_SPEC+channel+k;
							output[next_index] = 0;
						}
					#endif
					#if REPLACE_BROAD_ZERO == FALSE
						for(int l = START_CHANNEL; l < STOP_CHANNEL; l++){
							random = gasdev(&seed);
							next_index = (j+1)*CHANNELS_PER_SPEC+channel+k;
							output[next_index] = avg[next_index]+random*sd[next_index];
						}
					#endif

					#if PRINT_BROAD_STATS == TRUE
						printf("Output value: %d\n", output[current_index]);
						printf("Final window value: %f\n\n ", window[k]);
					#endif
				}
			}

			/*Recalculate stats for cleaned first thousand time samples (found in window array)*/
			for(int k = 1; k <= WINDOW_SIZE; k++){
				current_sample = (double)window[k-1];
				oldM = M;
				M = M + (current_sample-M)/k;
				S = S + (current_sample-M)*(current_sample-oldM);
			}
			V = S/(WINDOW_SIZE-1);
			SD = sqrt(V);
	
			/*Update statistics array*/
			/*
			for(int k = 0; k < WINDOW_SIZE; k++){
				current_index = k*channels_per_samp+channel;
				avg[current_index] = M;
				sd[current_index] = SD;
			}
			oldM = M;
			*/

			#if PRINT_BROAD_STATS == TRUE
				for(int k = 0; k < WINDOW_SIZE; k++){
					current_index = k*channels_per_samp+channel;
					printf("Current time sample:%d, and channel: %d\n", k, channel);
					printf("Input value: %f\n", input[current_index]);
					printf("Final mean value: %f\n", M);
					printf("Final standard deviation value: %f\n\n", S);
				}
			#endif
		}

		/*Begin calculation of rolling avg and sd*/
		if (j >= start_time + WINDOW_SIZE && j < num_samples - WINDOW_SIZE){
			#if PRINT_BROAD_STATS == TRUE
				printf("Current time sample:%ld\n", j);
			#endif

			/*Determine old, current, and next samples*/
			current_index = j*channels_per_samp+channel;
			next_index = (j+1)*channels_per_samp+channel;
			old_sample = window[0];
			current_sample = window[WINDOW_SIZE-1];
			next_sample = input[next_index];
			
			/*Update standard deviation and mean*/
			oldM = M;
			M = oldM + (current_sample-old_sample)/WINDOW_SIZE;
			V = V + (current_sample-old_sample)*(current_sample-M+old_sample-oldM)/(WINDOW_SIZE-1);
			SD = sqrt(V);

			/*Update statistics arrays*/
			avg[current_index] = M;
			sd[current_index] = SD;

   			/*Shift all elements in window by one*/ 
			memmove(&window[0], &window[1], (WINDOW_SIZE)*sizeof(double));
			window[WINDOW_SIZE-1] = next_sample;
			

			/*Set threshold limits*/
			thresh_top = M + thresh_limit*SD;
			thresh_bottom = M - thresh_limit*SD;

			#if PRINT_BROAD_STATS == TRUE
				printf("Mean: %f\n", M);
				printf("Standard deviation: %f\n", SD);
				printf("Input value: %f\n", input[current_index]);
				printf("Replaced input value (for stat calc): %f\n", current_sample);
				printf("Next value: %f\n", next_sample);
				printf("Threshold top: %f\n", thresh_top);
				printf("Threshold bottom: %f\n\n", thresh_bottom);
			#endif

			/*Replace bad values*/
			if(next_sample > thresh_top || next_sample < thresh_bottom){
				#if PRINT_BROAD_STATS1 == TRUE
					printf("******************************************************************\n");
					printf("HITTTT!!!!!!!!\n");
					printf("Current time sample: %ld\n", j);
					printf("*******************************************************************\n");
				#endif

				/*Increase counter for RFI detection*/
				counter++;
				
				/*Replace bad stat window value with noise*/	
				random = gasdev(&seed);
				/*while(random>1 || random<(-1)){
					random = gasdev(&seed);
				}*/
				window[WINDOW_SIZE-1] = M + SD*random;

				#if REPLACE_BROAD_ZERO == TRUE
					for(int k = START_CHANNEL; k < STOP_CHANNEL; k++){
						next_index = (j+1)*CHANNELS_PER_SPEC+channel+k;
						output[next_index] = 0;
					}
				#endif
				#if REPLACE_BROAD_ZERO == FALSE
					for(int k = START_CHANNEL; k < STOP_CHANNEL; k++){
						random = gasdev(&seed);
						next_index = (j+1)*CHANNELS_PER_SPEC+channel+k;
						output[next_index] = avg[next_index]+random*sd[next_index];
					}
				#endif

				#if PRINT_BROAD_STATS == TRUE
					//printf("Output value: %f\n\n", output[current_index]);
				#endif
			}
			/*
			for(int k = 0; k < WINDOW_SIZE; k++){
				printf("%f, ", window[k]);
			}
				printf("_________________________________\n");*/
		}
		
		/*For last WINDOW_SIZE samples, use stats calculated from the last WINDOW_SIZE samples before sample num_samples-WINDOW_SIZE*/
		if (j >= WINDOW_SIZE){
			current_index = j*channels_per_samp+channel;
			next_index = (j+1)*channels_per_samp+channel;
			next_sample = input[next_index];

			avg[current_index] = M;
			sd[current_index] = SD;

			if(next_sample > thresh_top || next_sample < thresh_bottom){
				#if PRINT_BROAD_STATS1 == TRUE
					printf("******************************************************************\n");
					printf("HITTTT!!!!!!!!\n");
					printf("Current time sample: %ld\n", j);
					printf("*******************************************************************\n");
				#endif

				/*Increase counter for RFI detection*/
				counter++;
				
				/*Replace bad stat window value with noise*/	
				random = gasdev(&seed);
				/*while(random>1 || random<(-1)){
					random = gasdev(&seed);
				}*/
				window[WINDOW_SIZE-1] = M + SD*random;

				#if REPLACE_BROAD_ZERO == TRUE
					for(int k = START_CHANNEL; k < STOP_CHANNEL; k++){
						next_index = (j+1)*CHANNELS_PER_SPEC+channel+k;
						output[next_index] = 0;
					}
				#endif
				#if REPLACE_BROAD_ZERO == FALSE
					for(int k = START_CHANNEL; k < STOP_CHANNEL; k++){
						random = gasdev(&seed);
						next_index = (j+1)*CHANNELS_PER_SPEC+channel+k;
						output[next_index] = avg[next_index]+random*sd[next_index];
					}
				#endif
			}

			#if PRINT_BROAD_STATS == TRUE
				printf("Mean: %f\n", M);
				printf("Standard deviation: %f\n", SD);
				printf("Input value: %f\n", input[current_index]);
				printf("Next value: %f\n", next_sample);
				printf("Threshold top: %f\n", thresh_top);
				printf("Threshold bottom: %f\n\n", thresh_bottom);
			#endif
		}
		
	}

	free(window);
	return counter;

}
