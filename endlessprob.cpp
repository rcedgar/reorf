#include "myutils.h"
#include "translator.h"

static byte GetRandLetter()
	{
	uint r = randu32();
	char c = "ACGT"[r%4];
	return (byte) c;
	}

static void GetRandSeq(byte *Seq, uint L)
	{
	for (uint i = 0; i < L; ++i)
		Seq[i] = GetRandLetter();
	}

static void EndlessProb(Translator &T, uint L, uint N,
  double &P1, double &P2)
	{
	byte *Seq = myalloc(byte, L);
	uint FoundCount1 = 0;
	uint FoundCount2 = 0;
	for (uint i = 0; i < N; ++i)
		{
		GetRandSeq(Seq, L);
		T.SetSeq(Seq, L, true);

		bool Plus, Minus;
		T.HasEndlessFrame(Plus, Minus);
		if (Plus || Minus)
			++FoundCount1;
		if (Plus && Minus)
			++FoundCount2;
		}
	myfree(Seq);

	P1 = double(FoundCount1)/N;
	P2 = double(FoundCount2)/N;
	}

static void Run1(uint Code, Translator &T, uint N, uint L_lo, uint L_hi,
  vector<uint> &Ls, vector<double> &P1s, vector<double> &P2s)
	{
	SetCode(Code);
	ProgressLog("Code %u\n", Code);

	Ls.clear();
	P1s.clear();
	P2s.clear();

	uint L_step = 50;
	uint L = L_lo;
	for (;;)
		{
		if (L == 500)
			L_step = 100;
		else if (L == 1000)
			L_step = 250;
		else if (L >= L_hi)
			break;

		for (uint i = 0; i < 3; ++i)
			{
			if ((L+i)%3 != 0)
				continue;
			double P1, P2;
			EndlessProb(T, L+i, N, P1, P2);
			Ls.push_back(L+i);
			P1s.push_back(P1);
			P2s.push_back(P2);
			Progress("Code %u L %u P1 %.3g P2 %.3g\n", Code, L+i, P1, P2);
			}

		L += L_step;
		}
	}

void LogP(double P)
	{
	Log("  %9.5f ", P);
	}

void cmd_endless_prob()
	{
	const string &SpecFileName = opt(endless_prob);

	uint L_lo = 150;
	uint L_hi = 2000;
	uint L_step = 50;

	uint N = 1000;
	if (optset_iters)
		N = opt(iters);

	Translator T;

	vector<uint> LSs;
	vector<uint> LMs;

	vector<double> PS1s;
	vector<double> PS2s;

	vector<double> PM1s;
	vector<double> PM2s;

	ProgressLog("N = %u (%s)\n", N, IntToStr(N));
	Run1(1, T, N, L_lo, L_hi, LSs, PS1s, PS2s);
	Run1(2, T, N, L_lo, L_hi, LMs, PM1s, PM2s);

	const uint M = SIZE(LSs);
	asserta(SIZE(LMs) == M);
	Log("\n");
	Log("\n");
	Log("           Standard code                        Mitochondrial code\n");
	Log("           =================================    ================================\n");
	Log("      L    P(Single)   P(Double)      P(D|S)    P(Single)  P(Double)      P(D|S)\n");
	//              123456789   123456789   123456789    123456789  123456789   123456789
	for (uint i = 0; i < M; ++i)
		{
		uint LS = LSs[i];
		uint LM = LMs[i];

		double PS1 = PS1s[i];
		double PM1 = PM1s[i];

		double PS2 = PS2s[i];
		double PM2 = PM2s[i];

		asserta(LS == LM);

		Log("   %4u  ", LS);
		LogP(PS1);
		LogP(PS2);
		if (PS1 > 0)
			LogP(PS2/PS1);
		else
			Log("  %9.9s ", ".");
		LogP(PM1);
		LogP(PM2);
		if (PM1 > 0)
			LogP(PM2/PM1);
		else
			Log("  %9.9s ", ".");
		Log("\n");
		}
	}
