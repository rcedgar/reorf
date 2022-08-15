#include "myutils.h"

byte TranslateCodon_Ints(uint i1, uint i2, uint i3);

static FILE *g_f;

static char GetLetter(uint i)
	{
	asserta(i < 4);
	return "ACGT"[i];
	}

static void CodonToStr(uint i1, uint i2, uint i3, string &Codon)
	{
	Codon.clear();
	Codon += GetLetter(i1);
	Codon += GetLetter(i2);
	Codon += GetLetter(i3);
	}

void DoCodonPair_OneDiff(uint i1, uint i2, uint i3,
  uint j1, uint j2, uint j3)
	{
	FILE *f = g_f;
	if (f == 0)
		return;

	int DiffCount = int(i1 != j1) + int(i2 != j2) + int(i3 != j3);
	asserta(DiffCount == 1);

	string Codon_i;
	string Codon_j;
	CodonToStr(i1, i2, i3, Codon_i);
	CodonToStr(j1, j2, j3, Codon_j);
	byte aa_i = TranslateCodon_Ints(i1, i2, i3);
	byte aa_j = TranslateCodon_Ints(j1, j2, j3);

	for (uint k = 0; k < 4; ++k)
		if (Codon_i[k] == Codon_j[k])
			{
			Codon_i[k] = tolower(Codon_i[k]);
			Codon_j[k] = tolower(Codon_j[k]);
			}

	fprintf(f, "%s", Codon_i.c_str());
	fprintf(f, "\t%s", Codon_j.c_str());
	fprintf(f, "\t%c", aa_i);
	fprintf(f, "\t%c", aa_j);
	if (aa_i == aa_j)
		fprintf(f, "\tSYN");
	else
		fprintf(f, "\tNON");
	fprintf(f, "\n");
	}

static void DoCodon(uint i1, uint i2, uint i3)
	{
	for (uint j = 1; j < 4; ++j)
		{
		DoCodonPair_OneDiff(i1, i2, i3, (i1+j)%4, i2, i3);
		DoCodonPair_OneDiff(i1, i2, i3, i1, (i2+j)%4, i3);
		DoCodonPair_OneDiff(i1, i2, i3, i1, i2, (i3+j)%4);
		}
	}

void cmd_mutable()
	{
	const string &OutputFileName = opt(mutable);
	g_f = CreateStdioFile(OutputFileName);

	for (uint i1 = 0; i1 < 4; ++i1)
		for (uint i2 = 0; i2 < 4; ++i2)
			for (uint i3 = 0; i3 < 4; ++i3)
				DoCodon(i1, i2, i3);

	CloseStdioFile(g_f);
	}
