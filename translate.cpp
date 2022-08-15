#include "myutils.h"
#include "alpha.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "translator.h"

void SeqToFasta(FILE *f, const char *Seq, unsigned L, const char *Label);

void GetFrames(vector<int> &Frames)
	{
	Frames.clear();
	if (optset_frame)
		{
		int Frame = atoi(opt(frame).c_str());
		Frames.push_back(Frame);
		}
	else
		{
		Frames.push_back(+1);
		Frames.push_back(+2);
		Frames.push_back(+3);
		Frames.push_back(-1);
		Frames.push_back(-2);
		Frames.push_back(-3);
		}
	}

static void WriteOrfsFasta(FILE *f, const Translator &T, const string &Label,
  uint &OrfIndex, bool Nucleo, bool EndlessOnly)
	{
	if (f == 0)
		return;

	string Sep = " ";
	if (optset_sep)
		Sep = opt(sep);

	const uint NtL = T.m_NtL;
	const vector<uint> &AaStarts = T.m_OrfAaStarts;
	const vector<uint> &AaLengths = T.m_OrfAaLengths;
	const vector<bool> &HasStops = T.m_OrfHasStops;
	const vector<bool> &HasMets = T.m_OrfHasMets;
	const vector<bool> &EndlessVec = T.m_OrfEndless;

	const uint N = SIZE(AaStarts);
	asserta(SIZE(AaLengths) == N);
	asserta(SIZE(HasStops) == N);

	for (uint i = 0; i < N; ++i)
		{
		int Frame = T.m_Frame;

		if (opt(plus_only) && Frame < 0)
			continue;
		if (opt(minus_only) && Frame > 0)
			continue;

		uint AaStart = AaStarts[i];
		uint AaLength = AaLengths[i];
		bool HasStop = HasStops[i];
		bool HasMet = HasMets[i];
		bool Endless = EndlessVec[i];

		if (EndlessOnly && !Endless)
			continue;

		string OrfSeq;
		if (Nucleo)
			T.GetNtSubStr(AaStart, AaLength, OrfSeq);
		else
			T.GetAaSubStr(AaStart, AaLength, OrfSeq);

		uint NtStart = T.GetNtPos(AaStart, 0);
		uint NtEnd = NtStart + 3*AaLength - 1;

		++OrfIndex;
		string OutLabel = Label;

		OutLabel += Sep;
		Psa(OutLabel, "orf%u", OrfIndex);

		OutLabel += Sep;
		Psa(OutLabel, "code%u", GetNCBICodeNr());

		OutLabel += Sep;
		Psa(OutLabel, "%u-%u/%u(%+d)",
		  NtStart + 1, NtEnd + 1, NtL, Frame);
		if (Endless)
			{
			if (NtL%3 == 0)
				OutLabel += "Ce";
			else
				OutLabel += "Co";
			if (!HasMet)
				OutLabel += "m";
			}

		if (Nucleo)
			SeqToFasta(f, OrfSeq.c_str(), 3*AaLength, OutLabel.c_str());
		else
			SeqToFasta(f, OrfSeq.c_str(), AaLength, OutLabel.c_str());
		}
	}

void RunTranslate(const string &InputFileName, bool GetOrfs)
	{
	const string &OutputFileName = opt(output);

	if (optset_aaout)
		Die("Use -output for aa seqs");

	FILE *fOut = CreateStdioFile(OutputFileName);
	FILE *fNtOut = CreateStdioFile(opt(ntout));

	FASTASeqSource SS;
	SS.Open(InputFileName);

	bool Circular = opt(circular);
	uint MinCodonCount = 30;
	if (optset_mincodons)
		MinCodonCount = opt(mincodons);

	bool SingleEndless = opt(single_endless);
	bool DoubleEndless = opt(double_endless);
	asserta(!(SingleEndless && DoubleEndless));
	if (SingleEndless || DoubleEndless)
		Circular = true;

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

		if (SingleEndless || DoubleEndless)
			{
			bool PlusHasEndless = false;
			bool MinusHasEndless = false;
			T.HasEndlessFrame(PlusHasEndless, MinusHasEndless);
			if (!PlusHasEndless && !MinusHasEndless)
				continue;
			if (DoubleEndless && (!PlusHasEndless || !MinusHasEndless))
				continue;
			if (SingleEndless && (PlusHasEndless && MinusHasEndless))
				continue;
			}

		for (uint k = 0; k < SIZE(Frames); ++k)
			{
			int Frame = Frames[k];
			if (!T.ValidFrame(Frame))
				continue;

			T.SetFrame(Frame);
			if (GetOrfs)
				{
				T.SetOrfs(MinCodonCount);
				WriteOrfsFasta(fOut, T, Label, OrfIndex,
				  false, SingleEndless || DoubleEndless);
				if (fNtOut != 0)
					WriteOrfsFasta(fNtOut, T, Label, OrfIndex,
					  true, SingleEndless || DoubleEndless);
				}
			else
				{
				string AaStr;
				T.GetAaStr(AaStr);

				string OutLabel;
				OutLabel = Label;
				Psa(OutLabel, " frame%+d", Frame);
				SeqToFasta(fOut, (const byte *) AaStr.c_str(), SIZE(AaStr),
				  OutLabel.c_str());
				}
			}
		}

	CloseStdioFile(fOut);
	CloseStdioFile(fNtOut);
	}

void cmd_translate()
	{
	RunTranslate(opt(translate), false);
	}

void cmd_getorfs()
	{
	RunTranslate(opt(getorfs), true);
	}
