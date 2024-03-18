#include "myutils.h"
#include "translator.h"

#define TRACE	0

static double Calcd(double p)
	{
	double x = (1 - 4.0*p/3.0);
	double d = -(3.0/4.0)*log(x);
	return d;
	}

void GetSynNonFracts(const char Codon[3], uint BaseNr,
  double &SynFract, double &NonFract)
	{
#if TRACE
	Log("\n");
	Log("  Codon %3.3s\n", Codon);
#endif
	char aa = TranslateCodon3(Codon);
	uint SynCount = 0;
	uint NonCount = 0;
	for (uint i = 0; i < 4; ++i)
		{
		char c = "ACGT"[i];
		if (Codon[BaseNr] == i)
			continue;

		char MutatedCodon[3];
		memcpy(MutatedCodon, Codon, 3);
		MutatedCodon[BaseNr] = c;
		char maa = TranslateCodon3(MutatedCodon);
		if (maa == '*')
			continue;
#if TRACE
		Log("   base %d/%c", BaseNr, c);
		Log(" mutated %3.3s", MutatedCodon);
		Log(" aa %c -> %c", aa, maa);
		Log("\n");
#endif
		if (maa == aa)
			++SynCount;
		else
			++NonCount;
		}
	SynFract = double(SynCount)/3;
	NonFract = double(NonCount)/3;

#if TRACE
	Log(" syn %.4f, nonsyn %.4f\n", SynFract, NonFract);
#endif
	}

static void DoPair(const string &QRowNt, const string &TRowNt,
  uint &SynCount, uint &NonCount,
  double &SumSynFract, double &SumNonFract)
	{
	Translator T;

	const uint ColCount = SIZE(QRowNt);
	asserta(ColCount%3 == 0);
	asserta(SIZE(TRowNt) == ColCount);
	const uint CodonCount = ColCount/3;

	if (opt(nogaps))
		{
		if (QRowNt.find('-') != string::npos)
			return;
		if (TRowNt.find('-') != string::npos)
			return;
		}

	string QRowAa;
	string TRowAa;

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
		if (qaa == '-' || taa == '-')
			continue;

		char QCodon[3];
		char TCodon[3];
		uint DiffCount = 0;
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			char qnt = QRowNt[3*CodonIndex+BaseNr];
			char tnt = TRowNt[3*CodonIndex+BaseNr];

			QCodon[BaseNr] = qnt;
			TCodon[BaseNr] = tnt;

			if (qnt != tnt)
				++DiffCount;
			}

		if (DiffCount > 1)
			continue;

		if (DiffCount == 1)
			{
			if (qaa == taa)
				++SynCount;
			else
				++NonCount;
			}

		double SumSynFract_ThisCodon = 0;
		double SumNonFract_ThisCodon = 0;
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			double SynFract, NonFract;
			GetSynNonFracts(QCodon, BaseNr, SynFract, NonFract);

			SumSynFract_ThisCodon += SynFract;
			SumNonFract_ThisCodon += NonFract;
			}
#if TRACE
		Log("  ns %.3g, nd %.3g\n", SumSynFract_ThisCodon, SumNonFract_ThisCodon);
#endif

		SumSynFract += SumSynFract_ThisCodon;
		SumNonFract += SumNonFract_ThisCodon;
		}
	}

void cmd_dsdn()
	{
	Die("TODO -- not working at all :-(");
#if 0
	const string TsvFN = opt(dsdn);

	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	vector<string> Fields;

	uint SynCount = 0;
	uint NonCount = 0;

	double SumSynFract = 0;
	double SumNonFract = 0;

	while (ReadLineStdioFile(f, Line))
		{
		if (Line[0] != 'A' || Line[1] != '\t')
			continue;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 5);
		const string &Query = Fields[1];
		const string &Target = Fields[2];
		const string &QRowNt = Fields[3];
		const string &TRowNt = Fields[4];

		DoPair(QRowNt, TRowNt, SynCount, NonCount,
		  SumSynFract, SumNonFract);

		DoPair(TRowNt, QRowNt, SynCount, NonCount,
		  SumSynFract, SumNonFract);
		}

	double pS = SynCount/SumSynFract;
	double pN = NonCount/SumNonFract;

	double dS = Calcd(pS);
	double dN = Calcd(pN);

	double dNdS = dN/dS;

	Log("SumSynFract %.1f, SumNonFract %.1f\n",
	  SumSynFract, SumNonFract);
	Log("SynCount %u, NonCount %u\n", SynCount, NonCount);

	Log("pN %.3g, pS %.3g\n", pN, pS);
	Log("dN %.3g, dS %.3g\n", dN, dS);
	Log("dN/dS %.3g\n", dNdS);

	CloseStdioFile(f);
#endif
	}
