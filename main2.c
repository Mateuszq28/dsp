#include <dsplib.h>


//Krok = (2*32768)/(48000/f) ≈ (f*22368)/16384 = (f*22368)>>14
//dla f=100Hz: Krok = STEP_SAW = round((2*32768)/(48000/f)) = round(136,5333) = 137
#define STEP_SAW 137
#define NUM_SAMPLES 5000
#define NUM_COEF 55
short samples[NUM_SAMPLES];
short samples2[NUM_SAMPLES];
//Bufor kołowy
//Wyzerowac i nie ruszać!!
int dbuffer[NUM_COEF+2];

//wspolczynniki filtru FIR
//Dlugosc filtru 55, liczba wspolczynnikow 55, rzad filtru 54
const short coefficients[] = {-1, -9, -18, -31, -48, -68,	-92, -116, -138, -153,
                -156, -142, -104, -39, 59, 190, 355, 551, 774, 1016,
                1268, 1520, 1760, 1975,	2156, 2293, 2378, 2407, 2378, 2293,
                2156, 1975, 1760, 1520, 1268, 1016, 774, 551, 355, 190,
                59, -39, -104, -142, -156, -153, -138, -116, -92, -68,
                -48, -31, -18, -9, -1};

//Filtr FIR
//input: wskaźnik do tablicy zawierającej próbki sygnału do przefiltrowania
//filter: wskaźnik do tablicy zawierającej współczynniki filtru
//output: wskaźnik do tablicy, w której zostaną zapisane wyniki filtracji
//numSamples: liczba próbek w tablicy
//numFilter: liczba współczynników filtru
void blockfir(short* input, const short* filter, short* output, int numSamples, int numFilter)
{
	int i;
	int j;
	long y;
	//Zaczynamy od ostatniej probki (najnowszej), zeby w przypadku podania
	//tej samej tablicy wejsciowej i wyjsciowej, móc ją nadpisać
	//(probek od konca nie bedziemy juz uzywac do obliczen)
	for(i = numSamples-1; i >= 0; i--)
	{
		for(j = 0; j < numFilter; j++)
		{
			if(i-j < 0) break;
			y = _smaci(y, filter[j], input[i-j]);
		}
		//Zaokrąglenie i zamiana z Q30 do Q15
		y = (y + (1<<14))>>15;
		output[i] = y;
	}
}

//Zwraca sumę wszystkich wartosci znajdujacych sie w tablicy buffer o dlugosci buflen
//Uwaga na przepelnienia!!!
int sumOfVector(int* buffer, unsigned int buflen)
{
	int i = 0;
	int sum = 0;
	for(i = 0; i < buflen; i++)
	{
		sum += buffer[i];
	}
	return sum;
}

//Funkcja wypelnia podany buffer o dlugosci buflen losowymi wartociami
//(szum biały)
void testfun(int* buffer, unsigned int buflen)
{
	rand16init();
	rand16((DATA*)buffer, buflen);
}

//Funkcja wypelnia podany buffer o dlugosci buflen sygnalem piloksztaltnym o stalym przyroscie step na kazda probke
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

//Funkcja wypelnia podany buffer o dlugosci buflen sygnalem prostokatnym
//sygnal jest tworzony na podstawie sygnalu piloksztaltnego o przyroscie step
//poziom wypelnienia jest sterowany parametrem threshold
//jest to prog decyzyjny - jesli sygnal piloksztaltny jest wiekszy od threshold
//syg prostokatny przybiera wartosci dodatnie, jesli mniejszy - ujemne
//dla poziomu wypelnienia 1/3 threshold = 10923 (wartosc w Q15)
//dla poziomu wypelnienia 50% threshold = 0
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

//Funkcja wypelnia podany buffer o dlugosci buflen sygnalem trojkatnym
//przyrost wynosi step*2 na każdą probke (jest on dodatni lub ujemny zaleznie od fazy)
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

//Funkcja wypelnia podany buffer o dlugosci buflen sygnalem sinusoidalnym
//faza sygnału odpowiada sygnalowi piloksztaltnemu o stalym przyroscie step na kazda probke
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

//Funkcja zeruje tablicę buffer o długosci buflen
void zeros(int* buffer, unsigned int buflen)
{
	int i;
	for(i = 0; i < buflen; i++)
		buffer[i] = 0;
}

void main(void)
{
	//GENERATORY
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
	//rand16init();
	//rand16((DATA*)samples, NUM_SAMPLES);

	//***************************************************
	//FILTRY FIR
	//int gain = sumOfVector((int*)coefficients, NUM_COEF); <- nie działa, bo przepelnienie

	//zad 2
	//testfun((int*)samples, NUM_SAMPLES);
	//blockfir(samples, coefficients, samples2, NUM_SAMPLES, NUM_COEF);
	//saw((int*)samples, NUM_SAMPLES, STEP_SAW);
	//blockfir(samples, coefficients, samples2, NUM_SAMPLES, NUM_COEF);

	//Testowanie czy można podac taki sam input i output -> można.
	//testfun((int*)samples, NUM_SAMPLES);
	//saw((int*)samples2, NUM_SAMPLES, STEP_SAW);
	//blockfir(samples, coefficients, samples, NUM_SAMPLES, NUM_COEF);
	//blockfir(samples2, coefficients, samples2, NUM_SAMPLES, NUM_COEF);

	//zad 3
	//testfun((int*)samples, NUM_SAMPLES);
	//zeros(dbuffer, NUM_COEF+2);
	//fir((DATA*)samples, (DATA*)coefficients, (DATA*)samples2, (DATA*)dbuffer, NUM_SAMPLES, NUM_COEF);
	//saw((int*)samples, NUM_SAMPLES, STEP_SAW);
	//zeros(dbuffer, NUM_COEF+2);
	//fir((DATA*)samples, (DATA*)coefficients, (DATA*)samples2, (DATA*)dbuffer, NUM_SAMPLES, NUM_COEF);

	//zad 4
	zeros(dbuffer, NUM_COEF+2);
	testfun((int*)samples, NUM_SAMPLES);
	int i;
	for(i = 0; i < NUM_SAMPLES; i++)
	{
		fir(&samples[i], (DATA*)coefficients, &samples2[i], (DATA*)dbuffer, 1, NUM_COEF);
	}
	saw((int*)samples, NUM_SAMPLES, STEP_SAW);
	for(i = 0; i < NUM_SAMPLES; i++)
	{
		fir(&samples[i], (DATA*)coefficients, &samples2[i], (DATA*)dbuffer, 1, NUM_COEF);
	}

	while (1); // do not exit
}
