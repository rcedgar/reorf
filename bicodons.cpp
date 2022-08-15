#include "myutils.h"
#include "translator.h"
#include "alpha.h"

/***
5-mer Mutation at position 3 only

123
GATTA
  321

4-mer Mutation at positions 2 or 3

123
GATT
 321
***/

static void CodeToStr(uint Code, uint k, string &Seq)
	{
	Seq.clear();
	for (uint i = 0; i < k; ++i)
		{
		uint b = (Code & 0x3);
		Seq += "ACGT"[b];
		Code >>= 2;
		}
	}

static void Domer(const string &Seq)
	{
	const uint k = SIZE(Seq);
	asserta(k == 4 || k == 5);

	char PlusCodon[3];
	char MinusCodon[3];

	PlusCodon[0] = Seq[0];
	PlusCodon[1] = Seq[1];
	PlusCodon[2] = Seq[2];

	if (k == 4)
		{
		MinusCodon[0] = g_CharToCompChar[Seq[3]];
		MinusCodon[1] = g_CharToCompChar[Seq[2]];
		MinusCodon[2] = g_CharToCompChar[Seq[1]];
		}
	else
		{
		MinusCodon[0] = g_CharToCompChar[Seq[4]];
		MinusCodon[1] = g_CharToCompChar[Seq[3]];
		MinusCodon[2] = g_CharToCompChar[Seq[2]];
		}

	bool InnerPal = (memcmp(PlusCodon, MinusCodon, 3) == 0);
	bool OuterPal = true;
	for (uint i = 0; i < k; ++i)
		{
		if (Seq[i] != g_CharToCompChar[Seq[k-i-1]])
			{
			OuterPal = false;
			break;
			}
		}

	char PlusAa = TranslateCodon3(PlusCodon);
	char MinusAa = TranslateCodon3(MinusCodon);

	uint LoMutPos = (k == 4 ? 1 : 2);
	for (uint MutPos = LoMutPos; MutPos <= 2; ++MutPos)
		{
		string MutSeq = Seq;
		char c = Seq[MutPos];
		for (uint i = 0; i < 4; ++i)
			{
			char Mutc = "ACGT"[i];
			if (Mutc == c)
				continue;

			MutSeq[MutPos] = Mutc;

			PlusCodon[0] = MutSeq[0];
			PlusCodon[1] = MutSeq[1];
			PlusCodon[2] = MutSeq[2];

			MinusCodon[0] = g_CharToCompChar[MutSeq[3]];
			MinusCodon[1] = g_CharToCompChar[MutSeq[2]];
			MinusCodon[2] = g_CharToCompChar[MutSeq[1]];

			char MutPlusAa = TranslateCodon3(PlusCodon);
			char MutMinusAa = TranslateCodon3(MinusCodon);

			Log("%s", Seq.c_str());
			Log("\t%u/%c/%c", MutPos, Seq[MutPos], MutSeq[MutPos]);
			Log("\t%s", MutSeq.c_str());
			Log("\t%c,%c", PlusAa, MinusAa);
			Log("\t%u%c", MutPos, Mutc);
			Log("\t%c,%c", MutPlusAa, MutMinusAa);

			if (PlusAa == '*' || MinusAa == '*' || MutPlusAa == '*' || MutMinusAa == '*')
				Log("\tstop");
			if (InnerPal)
				Log("\tinnerpal");
			if (OuterPal)
				Log("\touterpal");

			if (MutPlusAa == PlusAa && MutMinusAa == MinusAa)
				Log("\tsynboth");
			else if (MutPlusAa == PlusAa)
				Log("\tsynplus");
			else if (MutMinusAa == MinusAa)
				Log("\tsynminus");
			else
				Log("\tnonboth");
			Log("\n");
			}
		}
	Log("\n");
	}

void cmd_bicodons()
	{
	uint N4 = myipow(4, 4);
	uint N5 = myipow(4, 5);

	string Seq;
	for (uint Code = 0; Code < N4; ++Code)
		{
		CodeToStr(Code, 4, Seq);
		Domer(Seq);
		}

	for (uint Code = 0; Code < N5; ++Code)
		{
		CodeToStr(Code, 5, Seq);
		Domer(Seq);
		}
	}
