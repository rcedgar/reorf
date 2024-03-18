#include "myutils.h"
#include "translator.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"

void GetFrames(vector<int> &Frames);
void SeqToFasta(FILE *f, const char *Seq, unsigned L, const char *Label);

static void MatchMotif(const string &AaStr, const string &Motif, 
  uint &Pos, uint &Diffs)
	{
	Pos = 0;
	Diffs = UINT_MAX;
	const uint AL = SIZE(AaStr);
	const uint ML = SIZE(Motif);
	for (uint Pos2 = 0; Pos2 < AL; ++Pos2)
		{
		uint d = 0;
		for (uint k = 0; k < ML; ++k)
			{
			char a = AaStr[(Pos2 + k)%AL];
			char m = Motif[k];
			if (a != m)
				++d;
			}
		if (d < Diffs)
			{
			Diffs = d;
			Pos = Pos2;
			}
		}
	}

static void WriteNt(FILE *f, Translator &T, uint AaPos, const string &Label)
	{
	string Seq;
	T.GetNtStr(Seq);

	string RotatedSeq;
	const uint L = SIZE(Seq);
	for (uint i = 0; i < L; ++i)
		{
		char c = Seq[(3*AaPos + i)%L];
		RotatedSeq += c;
		}
	SeqToFasta(f, RotatedSeq.c_str(), L, Label.c_str());
	}

static void WriteAa(FILE *f, Translator &T, uint AaPos, const string &Label)
	{
	string Seq;
	T.GetAaStr(Seq);

	string RotatedSeq;
	const uint L = SIZE(Seq);
	for (uint i = 0; i < L; ++i)
		{
		char a = Seq[(AaPos + i)%L];
		RotatedSeq += a;
		}
	SeqToFasta(f, RotatedSeq.c_str(), L, Label.c_str());
	}

void cmd_rotate_motif()
	{
	const string &InputFileName = opt(rotate_motif);
	const string &Motif = opt(motif);

	FASTASeqSource SS;
	SS.Open(InputFileName);

	FILE *fAa = CreateStdioFile(opt(aaout));
	FILE *fNt = CreateStdioFile(opt(ntout));
	FILE *fTsv = CreateStdioFile(opt(tsvout));

	const bool Circular = true;
	Translator T;
	SeqInfo *SIQ = ObjMgr::GetSeqInfo();

	vector<int> Frames;
	GetFrames(Frames);
	for (;;)
		{
		bool Ok = SS.GetNext(SIQ);
		if (!Ok)
			break;

		const byte *Q = SIQ->m_Seq;
		const uint QL = SIQ->m_L;
		const string Label = string(SIQ->m_Label);

		T.SetSeq(Q, QL, Circular);
		uint OrfIndex = 0;

		uint BestDiffs = UINT_MAX;
		uint BestPos = UINT_MAX;
		int BestFrame = 0;
		for (uint k = 0; k < SIZE(Frames); ++k)
			{
			int Frame = Frames[k];
			if (!T.ValidFrame(Frame))
				continue;

			T.SetFrame(Frame);
			string AaStr;
			T.GetAaStr(AaStr);

			uint Pos, Diffs;
			MatchMotif(AaStr, Motif, Pos, Diffs);
			if (Diffs < BestDiffs)
				{
				BestDiffs = Diffs;
				BestPos = Pos;
				BestFrame = Frame;
				}
			}
		asserta(BestPos != UINT_MAX);
		asserta(BestFrame != 0);

		string NewLabel = Label + " rotated";
		T.SetFrame(BestFrame);
		WriteNt(fNt, T, BestPos, NewLabel);
		WriteAa(fAa, T, BestPos, NewLabel);

		Pf(fTsv, "%s	%+d	%u\n", Label.c_str(), BestFrame, BestPos);
		}

	CloseStdioFile(fAa);
	CloseStdioFile(fNt);
	CloseStdioFile(fTsv);
	}
