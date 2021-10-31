/* Convert LFP content from TDT format to dat format
 *
 * Usage:
 * tdt2dat <path_to_tank> <tank_name> <block_name> <var_name> <first_channel> <last_channel> <dat_filename>
 *
 * This is a very simple code to convert a set of channels stored in
 * TDT format to dat file.
 * At time I must implement some code more robust. For example:
 * Validate input variables;
 * Implement another code more flexible to, for example, convert only a sub set of the channels.
 *
 * Author: Nivaldo A P de Vasconcelos
 * Date: 2014May05*/


#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <errno.h>


/* TTank event header structure */
typedef struct tsqEventHeader
{
	uint64_t size_bytes;
	char type[3]; /* event type: snip, pdec, epoc etc */
	char version; /* event type: snip, pdec, epoc etc */
	char event_name[4]; /* event name: must be 4 chars, cast as a long */
	uint16_t channel; /* data acquisition channel */
	uint16_t num_channels; /* data acquisition channel */
	uint16_t sample_width_bytes;
	uint16_t reserved;
	uint8_t format; /* data format of event: byte, short, float (usually), double */
	uint8_t decimate;
	uint8_t rate; /* sampling frequency */
	uint64_t reserved2;
	uint32_t reserved3[2];
} tsqEventHeader_t;
#define NSAMPLES 1024*1024
//#define NSAMPLES 10
#define BYTES_PER_SAMPLE 4
//short int channel_map[64] = {1,8,2,7,3,58,6,4,5,9,16,10,15,11,14,12,13,24,17,18,23,19,22,20,21,25,32,26,27,31,30,28,29,33,41,34,39,35,38,36,37,40,49,43,46,42,47,45,44,48,57,51,54,50,55,53,52,56,64,59,63,62,61,60 };

// Channel map got from the physical mapping, by injecting signal in the adapter and checking in the TDT output signal.
//short int channel_map[64] = {1,8,2,7,3,6,4,5,9,16,10,15,11,14,12,13,17,24,18,23,19,22,20,21,25,32,26,31,27,30,28,29,33,41,34,39,35,38,36,37,40,49,43,46,42,47,45,44,48,57,51,54,50,55,53,52,56,64,59,63,58,62,61,60};
//short int channel_map[64] = {1,8,2,7,3,4,5,9,16,10,15,11,14,12,13,17,24,18,23,19,22,20,21,25,32,26,31,27,30,28,29,33,41,34,35,38,36,37,40,49,43,46,42,39,45,44,48,57,51,54,50,47,53,52,56,64,59,63,55,62,61,60,58,6};
// Old + channel map from Coherence Matrix 03Jun2014 
//short int channel_map[64] = {1,8,2,7,3,6,4,5,9,16,10,15,11,14,12,13,24,17,18,23,19,22,20,21,25,32,26,27,31,30,28,29,33,41,39,35,38,36,37,40,49,43,46,34,47,45,44,48,57,51,54,42,55,53,52,56,64,59,50,62,61,60,63,58};

short int channel_map[64]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64};

int main(int argc, char* argv[])
{

	FILE **sev = NULL, *dat = NULL;
	tsqEventHeader_t header;

	if (argc!=8) {
		fprintf (stderr,"\nWrong syntax:\n");
		fprintf (stderr,"Please, usage: tdt2dat <path_to_tanks> <tank_name> <block_name> <var_name> <first_channel> <last_channel> <dat_filename> \n");
		exit(0);
	}

	// Saving the input parameters
	char *path_to_tanks = argv[1];
	char *tank_name = argv[2];
	char *block_name = argv[3];
	char *var_name = argv[4];
	short int start_channel = atoi(argv[5]);
	short int last_channel = atoi(argv[6]);
	char *dat_filename = argv[7];


	short int tmp_int;
	if (last_channel<start_channel) {
		tmp_int = last_channel;
		last_channel = start_channel;
		start_channel = tmp_int;
	}

	short int num_channels = last_channel - start_channel + 1;


	char base_filename[1024-128];

	char sev_filename[1024];
	int i=0;
	float **buffer;
	int16_t  *sample_int;
	short int *channel;


	printf ("\nConverting TDT data file in the DAT format");
	printf ("\nTank path: %s",path_to_tanks);
	printf ("\nTank name: %s",tank_name);
	printf ("\nBlock name: %s",block_name);
	printf ("\nVariable name: %s",var_name);
	printf ("\nChannels: %d:%d",start_channel,last_channel);
	printf ("\nOutput DAT file: %s",dat_filename);
	printf ("\nProgress: ");



	channel = (short int *) malloc(sizeof(short int)*num_channels);
	for (i=0; i<num_channels; i++) {
		channel[i]=channel_map[(start_channel+i)-1];
		// channel[i]=start_channel+i; Original, before the channel mapping
	}

	// ---------------------------------------------------------------------------------
	// Opens the files: SEVs and DAT
	// ---------------------------------------------------------------------------------
		/* Allocate the memory to file pointers */
	sev = (FILE **) malloc(sizeof(FILE *)*num_channels);
	if (sev==NULL) {
		fprintf (stderr,"\nProblem to allocate memory. ErrorCode: SEV01 \n");
		return (1);
	}
		/* Open the SEV files */
	sprintf (base_filename,"%s/%s/%s/%s_%s_%s_Ch",path_to_tanks,tank_name,block_name,tank_name,block_name,var_name);
	for (i=0; i<num_channels; i++) {
		sprintf (sev_filename,"%s%d.sev",base_filename,channel[i]);
//		/printf ("\n%s",sev_filename);
		sev[i] = fopen(sev_filename, "r+b");
		if (sev[i]==NULL) {
			fprintf (stderr,"\nIt was not possible to open this file: %s\n",sev_filename);
			return (1);
		}
	}
		/* Open the DAT file */
	dat = fopen(dat_filename, "w+b");
	if (dat==NULL) {
		fprintf (stderr,"\nIt was not possible to open this file: %s; \t %s\n",dat_filename,strerror(errno));
		return (1);
	}
	// ---------------------------------------------------------------------------------


    // ---------------------------------------------------------------------------------
	// Allocate the buffers
	// ---------------------------------------------------------------------------------
		/* Allocate the pointers to array of buffers */
	buffer = (float **) malloc (sizeof(float *)*num_channels);
	if (buffer==NULL) {
		fprintf (stderr,"\nProblem to allocate memory. ErrorCode: SEV02 \n");
		return (1);
	}
		/* Allocate the array of buffer theirselves */
	for (i=0; i<num_channels; i++) {
		buffer[i] = (float *) malloc (sizeof(float)*NSAMPLES);
		if (buffer[i]==NULL) {
			fprintf (stderr,"\nProblem to allocate memory. ErrorCode: SEV03 \n");
			return (1);
		}
	}

	int count_bytes=0;
	int saved_bytes = 0;

	sample_int = (int16_t *) malloc (sizeof(int16_t)*(num_channels*NSAMPLES));
	if (sample_int==NULL) {
		fprintf (stderr,"\nProblem to allocate memory. ErrorCode: SEV04 \n");
		return (1);
	}
	// ---------------------------------------------------------------------------------

	// ---------------------------------------------------------------------------------
	// Reading operation
	// ---------------------------------------------------------------------------------
		/* Read the headers */
	for (i=0; i<num_channels; i++) {
		fread (&header,sizeof(tsqEventHeader_t),1,sev[i]);
	}
	//printf ("Firing rate: %2.4f",header.rate)*25000000/2^12/streamHeader.decimate)

		/* Read the data from channels */
	//int count = 0;
	int c=0;
	long int t=0;
//	float *tmp_data=NULL;
	long int total_saved_bytes = 0;
	int min_count_samples = INT_MAX;
	int num_samples = 0;
	while (!feof(sev[0])) {


			/* Fetches the samples for the next timestamp into buffer */
		saved_bytes = 0;
		min_count_samples = INT_MAX;
		for (i=0; i<num_channels; i++) {
			count_bytes = fread (buffer[i],BYTES_PER_SAMPLE,NSAMPLES,sev[i]);
			saved_bytes+= count_bytes;
			num_samples = count_bytes;
			if (num_samples<min_count_samples) {
				min_count_samples = num_samples;
			}
		}
			/* Convert to int16 */
		t = 0;
		for (c=0; c<min_count_samples; c++) {
			for (i=0; i<num_channels; i++) {
				sample_int[t++] = (int16_t) (*(buffer[i]+c)*1.0e6);
			}
		}
		fwrite(sample_int,2,t,dat);
		total_saved_bytes+=(saved_bytes*2);
		//printf ("%2.2f%% ",100.0*(total_saved_bytes/(float)(header.size_bytes*num_channels)));
		putchar ('.');
		fflush(stdout);
	}

	// ---------------------------------------------------------------------------------

	// ---------------------------------------------------------------------------------
	// Closing the files
	// ---------------------------------------------------------------------------------
	for (i=0; i<num_channels; i++) {
		fclose (sev[i]);
	}
	fclose (dat);
	printf ("\n[DONE]\n\n");
	// ---------------------------------------------------------------------------------


	return (0);

}


