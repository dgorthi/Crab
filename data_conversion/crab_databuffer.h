#include <stdint.h>
#include <stdio.h>

#define CACHE_ALIGNMENT         64
#define N_INPUT_BLOCKS          3 
#define N_OUTPUT_BLOCKS         3

#define N_CHAN_PER_PACK		512		//number of channels per packet
#define N_PACKETS_PER_SPEC	8		//number of packets per spectrum
#define N_BYTES_DATA_POINT	2		// number of bytes per datapoint in packet
#define N_POLS_CHAN		4		//number of polarizations per channel
#define N_BYTES_COUNTER		8		// number bytes of counter
#define N_CHANS_SPEC		N_CHAN_PER_PACK * N_PACKETS_PER_SPEC // not including the poles in spec. if we have 4 polarition, the total chans should time 4
#define DATA_SIZE_PACK		N_CHAN_PER_PACK * N_POLS_CHAN *  N_BYTES_DATA_POINT //doesn't include header for each packet size
#define PKTSIZE			N_CHAN_PER_PACK * N_POLS_CHAN *  N_BYTES_DATA_POINT + N_BYTES_COUNTER
#define N_BYTES_PER_SPEC	DATA_SIZE_PACK*N_PACKETS_PER_SPEC
#define N_SPEC_BUFF		128
#define BUFF_SIZE		N_SPEC_BUFF*N_BYTES_PER_SPEC
#define N_CHANS_BUFF		N_CHANS_SPEC*N_SPEC_BUFF*N_POLS_CHAN
