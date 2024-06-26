/*
 * Projekt z Zastosowania Procesorów Sygnałowych 2020
 * Projekt dla DSP TMS320C5535
 */


// Dołączenie wszelkich potrzebnych plików nagłówkowych
#include "usbstk5515.h"
#include "usbstk5515_led.h"
#include "aic3204.h"
#include "PLL.h"
#include "bargraph.h"
#include "oled.h"
#include "pushbuttons.h"
#include "dsplib.h"
#include "hamming.h"
#include "testsignal.h"


// Wzmocnienie wejścia w dB: 0 dla wejścia liniowego, 30 dla mikrofonowego
#define WZMOCNIENIE_WEJSCIA_dB 30

// Wybór sygnału do przetwarzania
// 0: próbki z wejścia mikrofonowego
// 1: biały szum
// 2: sygnał piłokształtny
// 3: sinus (stała częstotliwość)
// 4: sygnał z testsignal.h
#define SIGNAL 0

//liczba probek w sygnale wejsciowym (w przypadku generowania paczkami)
#define N 5000

//Krok = (2*32768)/(48000/f) ≈ (f*22368)/16384 = (f*22368)>>14
//dla f=100Hz: Krok = STEP_SAW = round((2*32768)/(48000/f)) = round(136,5333) = 137
#define STEP_SAW 137

#if SIGNAL != 0
	int input[N];
#endif

#define FILTER_BLOCK 512 //dlugosc bloku przetwarzania filterm w kroku 3
#define FFT_BLOCK 2048 //dlugosc bloku przetwarzania fft
#define NUM_COEF 55 //liczba wspolczynnikow filtru
//Bufor kołowy
//Wyzerowac i nie ruszać!!
int dbuffer[NUM_COEF+2];

//bufor pomocniczy na przetwarzane probki
//wybrac najwiekszy z FILTER_BLOCK, FFT_BLOCK
int samples[FFT_BLOCK];

//wspolczynniki filtru FIR
//Dlugosc filtru 55, liczba wspolczynnikow 55, rzad filtru 54
//okno Hamminga gain=32768
const short coefficients[] = {-1, -9, -18, -31, -48, -68,	-92, -116, -138, -153,
                -156, -142, -104, -39, 59, 190, 355, 551, 774, 1016,
                1268, 1520, 1760, 1975,	2156, 2293, 2378, 2407, 2378, 2293,
                2156, 1975, 1760, 1520, 1268, 1016, 774, 551, 355, 190,
                59, -39, -104, -142, -156, -153, -138, -116, -92, -68,
                -48, -31, -18, -9, -1};


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
	int j = 0;
	for(i = 0; i < limit; i = i+2)
	{
		//buffer[0] to usredniona wartosc sygnalu, skladowa stala
		re = _smpy(buffer[i], buffer[i]);
		im = _smpy(buffer[i+1], buffer[i+1]);
		buffer[j] = re + im; //mozemy dodac bez obawy przepelnienia, poniewaz sa to kwadraty liczb <1
		j++;
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

// Głowna procdura programu

void main(void) {

	int i; //uniwersalny licznik petli
	int freq; //pierwsze istotne maksimum widma (czestotliwosc)
	int index; //indeks pierwszego prazka w widmie
	int sig = 0; //licznik sygnalu wejsciowego

	// Inicjalizacja układu DSP
	USBSTK5515_init();			// BSL - układ uruchomieniowy
	pll_frequency_setup(100);	// ustawienie częstotliwości zegara DSP (PLL) na 100 MHz
	aic3204_hardware_init();	// kodek dźwięku AIC3204
	aic3204_init();				// jw.
	USBSTK5515_ULED_init();		// diody LED
	SAR_init_pushbuttons();		// przyciski
	oled_init();				// wyświelacz OLED 2x19 znaków

	// ustawienie częstotliwości próbkowania i wzmocnienia wejścia
	set_sampling_frequency_and_gain(48000L, WZMOCNIENIE_WEJSCIA_dB);

	// wypisanie komunikatu na wyświetlaczu
	// 2 linijki po 19 znaków, tylko wielkie angielskie litery
	oled_display_message("PROJEKT ZPS        ", "                   ");

	// 'krok' oznacza tryb pracy wybrany przyciskami
	unsigned int krok = 0;
	unsigned int poprzedni_krok = 9999;

	// zmienne do przechowywania wartości próbek
	Int16 lewy_wejscie;
	Int16 prawy_wejscie;
	Int16 lewy_wyjscie;
	Int16 prawy_wyjscie;
	Int16 mono_wejscie;

	// Przetwarzanie próbek sygnału w pętli
	while (1) {

		//sygnał wejsciowy
#if SIGNAL == 0
		// odczyt próbek audio, zamiana na mono
		aic3204_codec_read(&lewy_wejscie, &prawy_wejscie);
		mono_wejscie = (lewy_wejscie >> 1) + (prawy_wejscie >> 1);
#elif SIGNAL == 1
		if (sig == 0)
		{
			testfun(input, N);
		}
		mono_wejscie = input[sig];
		sig++;
		if (sig == N) sig = 0;
#elif SIGNAL == 2
		if (sig == 0)
		{
			saw(input, N, STEP_SAW);
		}
		mono_wejscie = input[sig];
		sig++;
		if (sig == N) sig = 0;

#elif SIGNAL == 3
		if (sig == 0)
		{
			saw(input, N, STEP_SAW);
			sine((DATA*)input, (DATA*)input, N);
		}
		mono_wejscie = input[sig];
		sig++;
		if (sig == N) sig = 0;

#elif SIGNAL == 4
		mono_wejscie = testsignal[sig];
		sig++;
		if (sig == N) sig = 0;
#endif

		// sprawdzamy czy wciśnięto przycisk
		// argument: maksymalna liczba obsługiwanych trybów
		krok = pushbuttons_read(4);
		if (krok == 0) // oba wciśnięte - wyjście
			break;
		else if (krok != poprzedni_krok) {
			// nastąpiła zmiana trybu - wciśnięto przycisk
			USBSTK5515_ULED_setall(0x0F); // zgaszenie wszystkich diód
			if (krok == 1) {
				// wypisanie informacji na wyświetlaczu
				oled_display_message("PROJEKT ZPS        ", "ORYGINALNY         ");
				// zapalenie diody nr 1
				USBSTK5515_ULED_on(0);
			} else if (krok == 2) {
				oled_display_message("PROJEKT ZPS        ", "FILTR              ");
				USBSTK5515_ULED_on(1);
				//zerowanie buforu kolowego
				zeros(dbuffer, NUM_COEF+2);
			} else if (krok == 3) {
				oled_display_message("PROJEKT ZPS        ", "FILTR BLOK         ");
				USBSTK5515_ULED_on(2);
				//zerowanie buforu kolowego i buforu pomocniczego we/wy
				zeros(dbuffer, NUM_COEF+2);
				zeros(samples, FILTER_BLOCK);
				//zerowanie licznika petli
				i = 0;
			} else if (krok == 4) {
				oled_display_message("PROJEKT ZPS        ", "WIDMO              ");
				USBSTK5515_ULED_on(3);
				//zerowanie licznika petli
				i = 0;
			}
			// zapisujemy nowo ustawiony tryb
			poprzedni_krok = krok;
		}


		// zadadnicze przetwarzanie w zależności od wybranego kroku

		if (krok == 1) {
			// tryb podstawowy - kopiowanie sygnału
			lewy_wyjscie = mono_wejscie;
			prawy_wyjscie = mono_wejscie;

		} else if (krok == 2) {
			//filtracja FIR
			fir((DATA*)mono_wejscie, (DATA*)coefficients, (DATA*)lewy_wyjscie, (DATA*)dbuffer, 1, NUM_COEF);
			prawy_wyjscie = lewy_wyjscie;

		} else if (krok == 3) {
			//filtracja blokowa FIR
			lewy_wyjscie = samples[i];
			prawy_wyjscie = lewy_wyjscie;
			samples[i] = mono_wejscie;
			if (i == FILTER_BLOCK-1)
			{
				fir((DATA*)samples, (DATA*)coefficients, (DATA*)samples, (DATA*)dbuffer, FILTER_BLOCK, NUM_COEF);
				i = 0;
			}
			else i++;

		} else if (krok == 4) {
			//analiza FFT
			samples[i] = mono_wejscie;
			lewy_wyjscie = mono_wejscie;
			prawy_wyjscie = lewy_wyjscie;
			if (i == FFT_BLOCK-1)
			{
				multiply(samples, (int*)hamming, FFT_BLOCK);
				rfft((DATA*)samples, FFT_BLOCK, SCALE);
				Amrfft(samples, FFT_BLOCK);
				index = maxIndex(samples, FFT_BLOCK, 987, 5);
				freq = freqIndex(index);
				i = 0;
			}
			else i++;

		}

		// zapisanie wartości na wyjście audio
		aic3204_codec_write(lewy_wyjscie, prawy_wyjscie);

	}

	// wciśnięcie obu przycisków jednocześnie kończy działanie pętli
	aic3204_disable();
	oled_display_message("KONIEC PRACY       ", "                   ");
	while(1); // nie wychodź z programu
}
