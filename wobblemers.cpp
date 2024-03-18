#include "myutils.h"
#include "fastaseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "translator.h"
#include "alpha.h"

const char *KmerCodeToStr(uint KmerCode, uint k, string &Str);
uint32 StrToKmerCode(const string &Seq);
uint32 SeqToKmerCode(const byte *Seq, uint k);
void TestKmerCode(uint k);
void LogDE(const char *Label, const byte *Seq, uint L);

static char GetPlusCodonFromWobblemer(const byte *Seq, char Codon[3])
	{
	Codon[0] = Seq[0];
	Codon[1] = Seq[1];
	Codon[2] = Seq[2];
	char aa = TranslateCodon3(Codon);
	return aa;
	}

static char GetMinusCodonFromWobblemer(const byte *Seq, char Codon[3])
	{
	Codon[0] = g_CharToCompChar[Seq[4]];
	Codon[1] = g_CharToCompChar[Seq[3]];
	Codon[2] = g_CharToCompChar[Seq[2]];
	char aa = TranslateCodon3(Codon);
	return aa;
	}

void WriteWobblemerCounts(FILE *f, const vector<uint> &Counts)
	{
	if (f == 0)
		return;

	const uint DictSize = 1024; // myipow(4, 5);
	asserta(SIZE(Counts) == DictSize);
	uint Total = 0;
	for (uint KmerCode = 0; KmerCode < DictSize; ++KmerCode)
		Total += Counts[KmerCode];

	double SumFreq = 0;
	for (uint KmerCode = 0; KmerCode < DictSize; ++KmerCode)
		{
		string Str;
		KmerCodeToStr(KmerCode, 5, Str);
		uint Count = Counts[KmerCode];
		double Freq = double(Count)/double(Total);
		SumFreq += Freq;
		if (Count > 0)
			{
			char PlusCodon[3];
			char MinusCodon[3];
			const byte *Wmer = (const byte *) Str.c_str();
			char aap = GetPlusCodonFromWobblemer(Wmer, PlusCodon);
			char aam = GetMinusCodonFromWobblemer(Wmer, MinusCodon);
			asserta(aap != '*');
			asserta(aam != '*');
			}
		fprintf(f, "%s\t%u\t%.5g\n", Str.c_str(), Count, Freq);
		}
	asserta(feq(SumFreq, 1));
	}

void CountWobblemers(SeqInfo *SI, vector<uint> &Counts)
	{
	const uint L = SI->m_L;
	Translator T;
	T.SetSeq(SI->m_Seq, L, true);
	int PlusEF = 0;
	int MinusEF = 0;
	T.GetEndlessFrames2(PlusEF, MinusEF);
	asserta(PlusEF > 0 && MinusEF < 0);
	uint Phase = GetPhase(PlusEF, MinusEF);
	if (Phase != 3)
		return;

	string PlusNtSeq;
	T.SetFrame(PlusEF);
	T.GetNtStr(PlusNtSeq);

	string MinusNtSeq;
	T.SetFrame(MinusEF);
	T.GetNtStr(MinusNtSeq);

	const byte *Q = SI->m_Seq;
	for (uint StartPos = uint(PlusEF)-1; StartPos + 5 < L; StartPos += 3)
		{
		const byte *Wmer = Q + StartPos;
		uint KmerCode = SeqToKmerCode(Wmer, 5);
		if (KmerCode != UINT32_MAX)
			{
			char PlusCodon[3];
			char MinusCodon[3];
			char aap = GetPlusCodonFromWobblemer(Wmer, PlusCodon);
			char aam = GetMinusCodonFromWobblemer(Wmer, MinusCodon);
			if (aap == '*')
				{
				LogDE(SI->m_Label, SI->m_Seq, SI->m_L);
				Log(">%s\n", SI->m_Label);
				Log("%*.*s\n", L, L, Q);
				Die("StartPos %u %5.5s has plus codon stop %3.3s",
				  StartPos, Wmer, PlusCodon);
				}
			if (aam == '*')
				{
				LogDE(SI->m_Label, SI->m_Seq, SI->m_L);
				Log(">%s\n", SI->m_Label);
				Log("%*.*s\n", L, L, Q);
				Die("StartPos %u %5.5s  has minus codon stop %3.3s",
				  StartPos, Wmer, MinusCodon);
				}

			asserta(KmerCode < SIZE(Counts));
			Counts[KmerCode] += 1;
			}
		}
	}

void cmd_wobblemers()
	{
	FASTASeqSource SS;
	SS.Open(opt(wobblemers));
	FILE *fOut = CreateStdioFile(opt(output));

	bool Circular = opt(circular);
	uint MinCodonCount = 30;
	if (optset_mincodons)
		MinCodonCount = opt(mincodons);

	bool SingleEndless = opt(single_endless);
	bool DoubleEndless = opt(double_endless);
	asserta(!(SingleEndless && DoubleEndless));
	if (SingleEndless || DoubleEndless)
		Circular = true;

	const uint DictSize = 1024; // myipow(4, 5);
	vector<uint> Counts(DictSize, 0);

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	for (;;)
		{
		bool Ok = SS.GetNext(SIQ);
		if (!Ok)
			break;

		CountWobblemers(SIQ, Counts);
		}

	WriteWobblemerCounts(fOut, Counts);
	}
