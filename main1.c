#include <dsplib.h>


//Krok = (2*32768)/(48000/f) ≈ (f*22368)/16384 = (f*22368)>>14
//Krok = STEP_SAW = round((2*32768)/(48000/f)) = round(136,5333) = 137
#define STEP_SAW 137
#define NUM_SAMPLES 5000
int samples[NUM_SAMPLES];
//int workBuf[NUM_SAMPLES];

void testfun(int* buffer, unsigned int buflen)
{
	rand16init();
	rand16((DATA*)buffer, buflen);
}

void saw(int* buffer, unsigned int buflen, int step)
{
	static int accumulator = 0;
	int i;
	for(i = 0; i < buflen; i++)
	{
		buffer[i] = accumulator;
		accumulator += step;
	}
}

void rect(int* buffer, unsigned int buflen, int step, int threshold)
{
	saw(buffer, buflen, step);
	int i;
	for(i = 0; i < buflen; i++)
	{
		if (buffer[i] < threshold)
			buffer[i] = -32768;
		else
			buffer[i] = 32767;
	}
}

void tri(int* buffer, unsigned int buflen, int step)
{
	saw(buffer, buflen, step);
	int i;
	int temp;
	for(i = 0; i < buflen; i++)
	{
		temp = buffer[i];
		temp = temp < 0 ? -temp : temp;
		buffer[i] = (temp - 16384)<<1;
	}
}

void sint(int* buffer, unsigned int buflen, int step)
{
	//mnożnik pi, wartoci <-1;1) - przemnożone przez pi daje faze
	saw(buffer, buflen, step);
	//stałe współczynniki ciągu Taylora przmnożone przez odpowiednie potęgi pi
	long a1 = 12868;
	long a3 = -21167;
	long a5 = 10445;
	long a7 = -2455;
	//argumenty do potegi 1, 2, 3, 5, 7
	int x1, x2, x3, x5, x7;
	long taylor;

	int i;
	int temp;
	int isPositive;
	for(i = 0; i < buflen; i++)
	{
		//faza trójkątna: 0 ... 0,5 ... 0 ... 0,5
		temp = buffer[i];
		if (temp >= 0) isPositive = 1;
		else isPositive = 0;
		temp = temp < 0 ? -temp : temp;
		if(temp > 16384) temp = 32767 - temp;

		//argumenty do potegi 1, 2, 3, 5, 7
		x1 = temp;
		x2 = _smpy(x1, x1);
		x3 = _smpy(x2, x1);
		x5 = _smpy(x3, x2);
		x7 = _smpy(x5, x2);

		//z szeregu Taylora
		taylor = a1*x1 + a3*x3 + a5*x5+ a7*x7;
		buffer[i] = (int)((taylor + (1<<11))>>12);
		if (isPositive == 0) buffer[i] = -buffer[i];
	}
}


void main(void)
{
	//zad 1
	//saw(samples, NUM_SAMPLES, STEP_SAW);

	//zad 2
	//rect(samples, NUM_SAMPLES, STEP_SAW, 0);
	//prog 1/3 -> (2*32768)*2/3-32768 = 10922 + 2/3
	//rect(samples, NUM_SAMPLES, STEP_SAW, 10923);

	//zad 3
	//tri(samples, NUM_SAMPLES, STEP_SAW);

	//zad 4
	//sint(samples, NUM_SAMPLES, STEP_SAW);

	//zad 5
	//saw(samples, NUM_SAMPLES, STEP_SAW);
	//sine((DATA*)samples, (DATA*)samples, NUM_SAMPLES);

	//zad 6
	rand16init();
	rand16((DATA*)samples, NUM_SAMPLES);

	while (1); // do not exit
}
