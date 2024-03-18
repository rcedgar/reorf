#include "myutils.h"
#include "translator.h"
#include "seqdb.h"

#define TRACE	0

static void DoPair(const string &QRowNt, const string &TRowNt,
  vector<uint> &Counts)
	{
	const uint ColCount = SIZE(QRowNt);
	asserta(ColCount%3 == 0);
	asserta(SIZE(TRowNt) == ColCount);
	const uint CodonCount = ColCount/3;

	string QRowAa;
	string TRowAa;

	Translator T;

	T.SetStr(QRowNt, false);
	T.SetFrame(1);
	T.GetAaStr(QRowAa);

	T.SetStr(TRowNt, false);
	T.SetFrame(1);
	T.GetAaStr(TRowAa);

	for (uint CodonIndex = 0; CodonIndex < CodonCount; ++CodonIndex)
		{
		char qaa = QRowAa[CodonIndex];
		char taa = TRowAa[CodonIndex];
		if (qaa == 'X' || taa == 'X')
			continue;
		if (qaa == '*' || taa == '*')
			continue;
		if (qaa == '-' || taa == '-')
			continue;

		if (CodonIndex+1 < CodonCount)
			{
			char q1 = QRowAa[CodonIndex+1];
			char t1 = TRowAa[CodonIndex+1];

			if (q1 == '-' || q1 == '*' || q1 == 'X')
				continue;
			if (t1 == '-' || t1 == '*' || t1 == 'X')
				continue;
			}

		uint DiffCount = 0;
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			char qnt = QRowNt[3*CodonIndex+BaseNr];
			char tnt = TRowNt[3*CodonIndex+BaseNr];

			if (qnt != tnt)
				++(Counts[BaseNr]);
			}
		}
	}

void cmd_freq3afa()
	{
	const string MSAFN = opt(freq3afa);

	SeqDB MSA;
	MSA.FromFasta(MSAFN, false);
	asserta(MSA.SetIsAligned());

	const uint SeqCount = MSA.GetSeqCount();

	vector<uint> Counts(3);
	for (uint i = 0; i < SeqCount; ++i)
		{
		string Rowi;
		MSA.GetSeqStr(i, Rowi);
		for (uint j = i + 1; j < SeqCount; ++j)
			{
			string Rowj;
			MSA.GetSeqStr(j, Rowj);

			DoPair(Rowi, Rowj, Counts);
			}
		}

	uint TotalCount = 0;
	for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
		{
		uint n = Counts[BaseNr];
		TotalCount += n;
		}

	for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
		{
		uint n = Counts[BaseNr];
		double Ratio = double(n)/TotalCount;
		ProgressLog("[%u]	%u	%.4f\n", BaseNr, n, Ratio);
		}

	double N12 = Counts[0] + Counts[1];
	double N3 = Counts[2];
	double Jiggle = 0.5*(N12 + 1)/(N3 + 1);
	ProgressLog("Jiggle %.3g\n", Jiggle);
	}
