#include "myutils.h"
#include "alpha.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "translator.h"

void SeqToFasta(FILE *f, const char *Seq, unsigned L, const char *Label);

uint GetPhase(int Endless_Plus, int Endless_Minus)
	{
#define x(p, m, Phase)	if (Endless_Plus == p && Endless_Minus == m) return Phase;

	x(+1, -1, 2)
	x(+1, -2, 3)
	x(+1, -3, 1)

	x(+2, -1, 3)
	x(+2, -2, 1)
	x(+2, -3, 2)

	x(+3, -1, 1)
	x(+3, -2, 2)
	x(+3, -3, 3)

#undef x

	asserta(false);
	return 0;
	}

static void GetFevStr(const SeqInfo *SIQ, string &Str)
	{
	Str.clear();

	const byte *Q = SIQ->m_Seq;
	const uint QL = SIQ->m_L;
	const string Label = string(SIQ->m_Label);

	Translator T;
	T.SetSeq(Q, QL, true);
	vector<int> Frames;
	vector<bool> HasMets;
	T.GetEndlessFrames(Frames, HasMets);
	const uint F = SIZE(Frames);
	if (F == 0)
		{
		int Frame = 0;
		uint OrfIndex = UINT_MAX;
		T.GetLongestOrf(Frame, OrfIndex);
		if (OrfIndex != UINT_MAX)
			{
			uint OrfLength = T.m_OrfAaLengths[OrfIndex];
			double Pct = (100.0*OrfLength)/double(QL);
			Ps(Str, "maxorf=%u/%.1f", OrfLength, Pct);
			}
		return;
		}

	uint NEpYesMet = 0;
	uint NEmYesMet = 0;
	uint NEpNoMet = 0;
	uint NEmNoMet = 0;
	for (uint i = 0; i < F; ++i)
		{
		int Frame = Frames[i];
		bool HasMet = HasMets[i];
		if (Frame > 0)
			{
			if (HasMet)
				++NEpYesMet;
			else
				++NEpNoMet;
			}
		else if (Frame < 0)
			{
			if (HasMet)
				++NEmYesMet;
			else
				++NEmNoMet;
			}
		else
			asserta(false);
		}

	Str = "maxorf=endless";
	if (NEpYesMet + NEmYesMet == 1)
		Str += "-single";
	else if (NEpYesMet + NEmYesMet == 2)
		Str += "-double";
	else if (NEpYesMet + NEmYesMet > 2)
		{
		Str += "-multi";
		uint MaxN = max(NEpYesMet, NEmYesMet);
		uint MinN = min(NEpYesMet, NEmYesMet);
		Psa(Str, "/%u/%u", MaxN, MinN);
		}

	if (NEpNoMet + NEmNoMet > 0)
		{
		uint MaxN = max(NEpNoMet, NEmNoMet);
		uint MinN = min(NEpNoMet, NEmNoMet);
		Psa(Str, "-startless/%u/%u", MaxN, MinN);
		}
	}

void cmd_getorfs_serratus()
	{
	const string &InputFileName = opt(getorfs_serratus);
	const string &OutputFileName = opt(fevout);

	if (optset_output)
		Die("Use -fevout");

	FILE *fOut = CreateStdioFile(OutputFileName);

	FASTASeqSource SS;
	SS.Open(InputFileName);
	uint64 FileSize = GetStdioFileSize64(SS.m_LR.m_f);

	uint MinCodonCount = 30;
	if (optset_mincodons)
		MinCodonCount = opt(mincodons);

	bool Circular = true;

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	asserta(!optset_frame);

	uint SeqCount = 0;
	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SIQ);
		if (!Ok)
			break;

		++SeqCount;
		if (SeqCount%1000 == 0)
			{
			uint64 FilePos = GetStdioFilePos64(SS.m_LR.m_f);
			double Fract = double(FilePos)/double(FileSize);
			uint Mills = uint(Fract*1000);
			if (Mills == 0)
				Mills = 1;
			if (Mills >= 999)
				Mills = 998;
			ProgressStep(Mills, 1000, "Processing (%s seqs)", IntToStr(SeqCount));
			}

		string FevStr;
		GetFevStr(SIQ, FevStr);

		if (FevStr != "" && fOut != 0)
			fprintf(fOut, "%s\t%s\n", SIQ->m_Label, FevStr.c_str());
		}
	ProgressStep(999, 1000, "Processing (%s seqs)", IntToStr(SeqCount));

	CloseStdioFile(fOut);
	}
