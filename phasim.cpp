#include "myutils.h"
#include "translator.h"

uint GetPhase(int Endless_Plus, int Endless_Minus);
byte TranslateCodon(byte b1, byte b2, byte b3);

byte GetRandChar()
	{
	return "ACGT"[randu32()%4];
	}

static void GetRandNonStopCodon(byte *Codon)
	{
	for (uint k = 0; k < 1000; ++k)
		{
		byte b1 = GetRandChar();
		byte b2 = GetRandChar();
		byte b3 = GetRandChar();
		byte aa = TranslateCodon(b1, b2, b3);
		if (aa != '*')
			{
			Codon[0] = b1;
			Codon[1] = b2;
			Codon[2] = b3;
			return;
			}
		}
	asserta(false);
	}

void GetRandEndlessOrf(byte *Seq, uint CodonCount)
	{
	for (uint i = 0; i < CodonCount; ++i)
		GetRandNonStopCodon(Seq + i*3);
	}

void cmd_phasim()
	{
	FILE *fOut = CreateStdioFile(opt(phasim));
	const uint IterCount = 10000;

	Translator T;

	const uint MinCodonCount = 150;
	const uint MaxCodonCount = 500;
	const uint d = MaxCodonCount - MinCodonCount;
	const bool Circular = true;
	byte *Seq = myalloc(byte, MaxCodonCount*3);
	vector<uint> StopCounts(3, 0);
	vector<uint> LetterCounts(3, 0);
	for (uint Iter = 0; Iter < IterCount; ++Iter)
		{
		ProgressStep(Iter, IterCount, "Simulating");
		uint CodonCount = MinCodonCount + randu32()%d;
		asserta(CodonCount >= MinCodonCount);
		asserta(CodonCount <= MaxCodonCount);
		GetRandEndlessOrf(Seq, CodonCount);

		T.SetSeq(Seq, CodonCount*3, Circular);
		for (uint uFrame = 1; uFrame <= 3; ++uFrame)
			{
			T.SetFrame(-int(uFrame));
			string AaStr;
			T.GetAaStr(AaStr);
			const uint AL = SIZE(AaStr);
			LetterCounts[uFrame-1] += AL;
			for (uint i = 0; i < AL; ++i)
				{
				if (AaStr[i] == '*')
					StopCounts[uFrame-1] += 1;
				}
			}
		}

	for (uint uFrame = 1; uFrame <= 3; ++uFrame)
		{
		uint n = StopCounts[uFrame-1];
		uint M = LetterCounts[uFrame-1];
		int iFrame = -int(uFrame);
		uint Phase = GetPhase(1, iFrame);
		double Freq = double(n)/M;

		ProgressLog("Frame -%u, phase %u, stops %u / %u, freq %.3g\n",
		  uFrame, Phase, n, M, Freq);
		Pf(fOut, "Frame -%u, phase %u, stops %u / %u, freq %.3g\n",
		  uFrame, Phase, n, M, Freq);
		}
	}
