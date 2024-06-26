#include <dsplib.h>
#include <stdio.h>
#include "testsignal.h"
#include "hamming.h"


//Krok = (2*32768)/(48000/f) ≈ (f*22368)/16384 = (f*22368)>>14
//dla f=100Hz: Krok = STEP_SAW = round((2*32768)/(48000/f)) = round(136,5333) = 137
#define STEP_SAW 137
#define NUM_SAMPLES 5000
#define NUM_COEF 55

//DO USTAWIENIA FFT
#define N 2048 //dlugosc okna przetwarzania
#pragma DATA_SECTION (buffer_fft, ".input")
//int buffer_fft[N];

int samples[N];
int samples2[N];
//Bufor kołowy
//Wyzerowac i nie ruszać!!
int dbuffer[NUM_COEF+2];


//wspolczynniki filtru FIR
//Dlugosc filtru 55, liczba wspolczynnikow 55, rzad filtru 54
//okno Hamminga gain=32768
const short coefficients[] = {-1, -9, -18, -31, -48, -68,	-92, -116, -138, -153,
                -156, -142, -104, -39, 59, 190, 355, 551, 774, 1016,
                1268, 1520, 1760, 1975,	2156, 2293, 2378, 2407, 2378, 2293,
                2156, 1975, 1760, 1520, 1268, 1016, 774, 551, 355, 190,
                59, -39, -104, -142, -156, -153, -138, -116, -92, -68,
                -48, -31, -18, -9, -1};
//okno Blackmana gain=32768
const short coefficients2[] = {0, 0, -1, -3, -8, -15, -26, -38, -53, -66,
				-75, -75, -60, -24, 39, 136, 269, 441, 650, 891,
				1154, 1429, 1699, 1950, 2166, 2332, 2436, 2472, 2436, 2332,
				2166, 1950, 1699, 1429, 1154, 891, 650, 441, 269, 136,
				39, -24, -60, -75, -75, -66, -53, -38, -26, -15,
				-8, -3, -1, 0, 0};

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

//Odwraca kolejnosc probek w tablicy buffer o dlugosci buflen
void inverse(int* buffer, unsigned int buflen)
{
	int limit = buflen>>1;
	int i;
	int temp;
	for(i = 0; i < limit; i++)
	{
		temp = buffer[i];
		buffer[i] = buffer[buflen-1-i];
		buffer[buflen-1-i] = temp;
	}
}

//Kopiuje buflen pierwszych probek z wektora input do wektora output
void copy(int* input, int* output, unsigned int buflen)
{
	int i;
	for(i = 0; i < buflen; i++)
		output[i] = input[i];
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

//Oblicza widmo amplitudowe na podstawie zespolonego widma
//Nadpisuje probki widma z buffer
//buflen musi byc potegą liczby 2
//STRUKTURA BUFFER:
//dwie pierwsze probki to: skladowa stala, skladowa Nyquista: y(0)Re, y(nx/2)im
//pozostałe: y(1)Re, y(1)Im, y(2)Re, y(2)Im...
void Amrfft(int* buffer, unsigned int buflen)
{
	//uzyc w przypadku liczenia dokladnej amplitudy na probke
	//opcja 2. w skalowaniu
	//int pow;
	//log_2((DATA*)buflen, (LDATA*)pow, 1);

	int limit = buflen>>1;
	int re, im;
	int i;
	for(i = 0; i < limit; i++)
	{
		//buffer[0] to usredniona wartosc sygnalu, skladowa stala
		re = _smpy(buffer[2*i], buffer[2*i]);
		im = _smpy(buffer[2*i+1], buffer[2*i+1]);
		buffer[i] = re + im; //mozemy dodac bez obawy przepelnienia, poniewaz sa to kwadraty liczb <1
	}

	sqrt_16((DATA*)buffer, (DATA*)buffer, limit); //po tej operacji w temp mamy modul widma (jena polowa)

	//Skalowanie
	for(i = 0; i < limit; i++)
	{
		//1. musimy uwzglednic drugą (symetryczną) czesc widma
		buffer[i] = buffer[i]<<1;

		//2. ponizsze jesli chcielibysmy znac dokladna amplitude przypadającą na próbke
		//tutaj za mala precyzja i ucina 6 prazek z widma (dla danych klarnet)
		//buffer[i] = buffer[i]>>(pow-1);
	}

	//zerwowanie nadmiarowych probek
	for(i = limit; i < buflen; i++)
		buffer[i] = 0;
}

//Mnozy bufor inOut przez in i zapisuje w inOut
//Oba bufory mają dlugosc buflen
void multiply(int* inOut, int* in, unsigned int buflen)
{
	int i;
	for(i = 0; i < buflen; i++)
		inOut[i] = _smpy(inOut[i], in[i]);
}

//Zwraca indeks pierwszego maksimum lokalnego odrozniajacego sie od szumu z tablicy buffer o dlugosci buflen
//maksimum lokalne musi byc wieksze od threshold (odrzucanie szumu)
//maksymalna szerokosc prążka - maxWidth
int maxIndex(int* buffer, unsigned int buflen, int threshold, int maxWidth)
{
	int i;
	int deltaO = 0; //poprzednia delta
	int deltaN = 0; //obecna delta
	int summit = 0; //licznik szerokosci szczytu
	for(i = 2; i < buflen; i++)
	{
		deltaO = buffer[i-1] - buffer[i-2];
		deltaN = buffer[i] - buffer[i-1];

		if(deltaO > 0)
		{
			if(summit == 0)
			{
				if(deltaN == 0) summit++;
				if (deltaN < 0)
				{
					if(buffer[i-1] > threshold) return i-1;
				}
			}
		}

		if(deltaO == 0)
		{
			if(deltaN > 0 && summit > 0) summit = 0;
			if(summit > 0)
			{
				if(deltaN == 0) summit++;
				if(deltaN < 0)
				{
					if(summit <= maxWidth && buffer[i-1] > threshold) return i-1-(summit>>1);
					else summit = 0;
				}
			}
		}
	}
	return 0;
}


//Zwraca czestotliwosc prazka widma o indeksie index
//Do obliczen zalozono:
//szybkosc probkowania 48 kHz
//staly rozmiar FFT 2048
int freqIndex(int index)
{
	if (index <= 87) return (index*375)>>4;
	else if (index <= 174) return (((index*125)>>1)*3)>>3;
	else if (index <= 436) return (((index*75)>>3)*5)>>1;
	else if (index <= 699) return ((((index*25)>>1)*3)>>2)*5>>1;
	else if (index <= 1048) return ((((index*25)>>2)*5)>>2)*3;
	else if (index <= 1310) return ((((index*25)>>2)*3)>>2)*5;
	else if (index <= 1398) return ((((index*15)>>2)*5)>>2)*5;
	else return ((((index*15)>>3)*5)>>1)*5;
}

//Drukuje czestotliwosc prazka o indksie indeks
//Do obliczen zalozono:
//szybkosc probkowania 48 kHz
//staly rozmiar FFT 2048
void printFreq(int index)
{
	int freq = freqIndex(index);
	if (freq < 0)
	{
		freq -= 32767;
		printf("Czestotliwosc prazka o indeksie %d wynosi: 32767+%d\n", index, freq);
		//printf("=32767+%d\n", freq); //do porownania w excelu
	}
	else
	{
		printf("Czestotliwosc prazka o indeksie %d wynosi: %d\n", index, freq);
		//printf("%d\n", freq); //do porownania w excelu
	}
}

void main(void)
{
	//GENERATORY
	//zad 1 Generowanie sygnału piłokształtnego
	//saw(samples, NUM_SAMPLES, STEP_SAW);

	//zad 2 Generowanie sygnału prostokątnego
	//rect(samples, NUM_SAMPLES, STEP_SAW, 0);
	//prog 1/3 -> (2*32768)*2/3-32768 = 10922 + 2/3
	//rect(samples, NUM_SAMPLES, STEP_SAW, 10923);

	//zad 3 Generowanie sygnału trójkątnego
	//tri(samples, NUM_SAMPLES, STEP_SAW);

	//zad 4 Generowanie sygnału sinus za pomocą szeregu Taylora
	//sint(samples, NUM_SAMPLES, STEP_SAW);

	//zad 5 Generowanie sygnału sinus za pomocą funkcji z DSPLIB
	//saw(samples, NUM_SAMPLES, STEP_SAW);
	//sine((DATA*)samples, (DATA*)samples, NUM_SAMPLES);

	//zad 6 Generowanie białego szumu
	//rand16init();
	//rand16((DATA*)samples, NUM_SAMPLES);

	//***************************************************
	//FILTRY FIR
	//int gain = sumOfVector((int*)coefficients, NUM_COEF); <- nie działa, bo przepelnienie

	//zad 2 Własna implementacja filtru FIR
	//testfun((int*)samples, NUM_SAMPLES);
	//blockfir(samples, coefficients, samples2, NUM_SAMPLES, NUM_COEF);
	//saw((int*)samples, NUM_SAMPLES, STEP_SAW);
	//blockfir(samples, coefficients, samples2, NUM_SAMPLES, NUM_COEF);

	//Testowanie czy można podac taki sam input i output -> można.
	//testfun((int*)samples, NUM_SAMPLES);
	//saw((int*)samples2, NUM_SAMPLES, STEP_SAW);
	//blockfir(samples, coefficients, samples, NUM_SAMPLES, NUM_COEF);
	//blockfir(samples2, coefficients, samples2, NUM_SAMPLES, NUM_COEF);

	//zad 3 Filtracja blokowa sygnału za pomocą funkcji z DSPLIB
	//testfun((int*)samples, NUM_SAMPLES);
	//zeros(dbuffer, NUM_COEF+2);
	//fir((DATA*)samples, (DATA*)coefficients, (DATA*)samples2, (DATA*)dbuffer, NUM_SAMPLES, NUM_COEF);
	//saw((int*)samples, NUM_SAMPLES, STEP_SAW);
	//zeros(dbuffer, NUM_COEF+2);
	//fir((DATA*)samples, (DATA*)coefficients, (DATA*)samples2, (DATA*)dbuffer, NUM_SAMPLES, NUM_COEF);

	//zad 4 Filtracja sygnału próbka po próbce za pomocą funkcji z DSPLIB
	//zeros(dbuffer, NUM_COEF+2);
	//testfun((int*)samples, NUM_SAMPLES);
	//int i;
	//for(i = 0; i < NUM_SAMPLES; i++)
	//{
	//	fir(&samples[i], (DATA*)coefficients, &samples2[i], (DATA*)dbuffer, 1, NUM_COEF);
	//}
	//saw((int*)samples, NUM_SAMPLES, STEP_SAW);
	//for(i = 0; i < NUM_SAMPLES; i++)
	//{
	//	fir(&samples[i], (DATA*)coefficients, &samples2[i], (DATA*)dbuffer, 1, NUM_COEF);
	//}

	//TESTOWANIE FILTRU FIR NA ODWRÓCONYM SYGNALE
	//int x = testsignal[0]; //aby prążek byl widoczny w "Expressions"
	//blockfir((short*)testsignal, coefficients, samples, 2048, NUM_COEF);
	//copy((int*)testsignal, (int*)samples, 2048);
	//inverse((int*)samples, 2048);
	//blockfir(samples, coefficients, samples, 2048, NUM_COEF);
	//inverse((int*)samples, 2048);

	//***************************************************
	//WIDMO

	//zad 1 Przygotowanie sygnału testowego
	int x = testsignal[0]; //aby prążek byl widoczny w "Expressions"

	//zad 2 Analiza widmowa sygnału
	//copy((int*)testsignal, buffer_fft, N);
	//multiply(buffer_fft, (int*)hamming, N);
	//rfft((DATA*)buffer_fft, N, SCALE);
	//Amrfft(buffer_fft, N);

	//zad 3 Zastosowanie funkcji okna
	x = hamming[0];

	//zad 4 Wyszukiwanie maksimów w widmie
	//copy((int*)testsignal, buffer_fft, N);
	//rfft((DATA*)buffer_fft, N, SCALE);
	//Amrfft(buffer_fft, N);

	//Wypisanie pierwszego prążka
	//int index = maxIndex(buffer_fft, N<<2, 987, 5);
	//printf("Pierwszy prazek ma indeks: %d \n", index);
	//int freq = freqIndex(index);
	//printFreq(index);

	//Sprawdzenie pozostalych czestotliwosci
	//int i;
	//printf("\n");
	//for(i = 0; i < 2048; i++)
	//	printFreq(i);

	//zad 5 Znajdywanie częstotliwości podstawowej metodą autokorelacji
	int i;
	for (i = 0; i < N; i++)
		samples[i] = testsignal[i]>>4;
	acorr((DATA*)samples, (DATA*)samples2, 2048, 2048, bias);
	int index = maxIndex(samples2, N, 0, 5);
	printf("Pierwsze maksimum koleracji wystepuje dla indeksu: %d \n", index);


	while (1); // do not exit
}
