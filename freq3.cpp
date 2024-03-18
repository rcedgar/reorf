#include "myutils.h"
#include "translator.h"

#define TRACE	0

static void DoPair(const string &QRowNt, const string &TRowNt,
  vector<uint> &Counts)
	{
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
		asserta(qaa != '*' && taa != '*');
		if (qaa == '-' || taa == '-')
			continue;

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

void cmd_freq3()
	{
	const string TsvFN = opt(freq3);

	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	vector<string> Fields;

	vector<uint> Counts(3);
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

		DoPair(QRowNt, TRowNt, Counts);
		}

	CloseStdioFile(f);

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
