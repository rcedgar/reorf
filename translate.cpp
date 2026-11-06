#include "myutils.h"
#include "alpha.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "translator.h"

// -getorfs and -getorfs_serratus
//
// Both commands find ORFs and can write the same files. Defaults differ so
// that existing command lines keep their results.
//
//   -getorfs            linear, unless -circular, -single_endless, or
//                       -double_endless is set
//   -getorfs_serratus   circular, unless -linear is set
//
// Outputs, each written only when its filename is set. Sequence files list
// every qualifying ORF unless -longest_orf is set.
//   -output    amino-acid ORF FASTA. The terminal stop codon is omitted.
//   -ntout     nucleotide ORF FASTA, same records and labels as -output.
//   -fevout    Serratus TSV: input label, tab, maxorf=...
//              Without -longest_orf, stop-free frames are classified
//              (maxorf=endless-single, -double, -multi/..., -startless/...).
//              maxorf=LEN/PCT is written only when no selected frame is
//              endless. LEN is the single longest ORF across those frames,
//              not one ORF per strand. PCT is 100 * amino-acid length /
//              nucleotide length. With -longest_orf the line is always
//              maxorf=LEN/PCT for the one selected ORF, including when a
//              frame is endless. A sequence with nothing to report produces
//              no line.
//   -alnout    Forward nucleotide sequence on one row, then one row per
//              reported ORF.
//              Each amino acid is the single-letter code placed on the middle
//              base of its codon ('*' = stop). The terminal stop is shown
//              here even though FASTA omits it. Rows are split into blocks
//              of 100 nucleotides unless -fasta_cols is set (0 = one block).
//              FASTA line length stays 80 when -fasta_cols is omitted.
//              Two spaces and a label follow each row: the input FASTA label
//              on the nucleotide row, the FASTA ORF label on an ORF row.
//              A blank line separates blocks.
//
// Other options:
//   -frame              report this frame only (+1,+2,+3,-1,-2,-3)
//   -plus_only          ignore minus frames
//   -minus_only         ignore plus frames
//   -longest_orf        report only the longest ORF in every output file.
//                       Greater codon length wins. Equal length prefers the
//                       plus strand. The same strand keeps the earlier ORF
//                       (+1,+2,+3,-1,-2,-3, then the order found in the frame).
//                       orfN stays the index from the full list.
//   -mincodons          minimum ORF length in codons (default 30)
//   -sep                label field separator (default space)
//   -circular           force circular topology
//   -linear             force linear topology
//   -single_endless     keep sequences with an endless ORF on exactly one
//                       strand; write only endless ORFs; forces circular
//   -double_endless     keep sequences with an endless ORF on both strands;
//                       write only endless ORFs; forces circular
//   -code               NCBI genetic code (default 1)
//
// -circular and -linear together is an error. -linear with -single_endless
// or -double_endless is an error. Endless-sequence filtering inspects every
// valid frame; ORF records are then taken from the selected frames.
//
// ORF label:
//   <input><sep>orfN<sep>codeC<sep>start-end/len(+frame)[Ce|Co][m]
// Ce = circular endless, length divisible by 3
// Co = circular endless, length not divisible by 3
// m  = endless ORF with no methionine

void SeqToFasta(FILE *f, const char *Seq, unsigned L, const char *Label);
void GetFevStr(const SeqInfo *SIQ, bool Circular, const vector<int> &Frames,
  uint MinCodonCount, bool PlusOnly, bool MinusOnly, string &Str);

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

struct OrfHit
	{
	string Label;
	string AaSeq;
	string NtSeq;
	string AlnRow;
	uint AaLength = 0;
	int Frame = 0;
	};

static string MakeAlnRow(const Translator &T, uint AaStart, uint AaLength,
  bool HasStop)
	{
	const uint NtL = T.m_NtL;
	string Row(NtL, ' ');
	for (uint Off = 0; Off < AaLength; ++Off)
		{
		uint Pos = T.GetNtPos(AaStart + Off, 1);
		asserta(Pos < NtL);
		asserta(Row[Pos] == ' ');
		Row[Pos] = (char) T.GetAa(AaStart + Off);
		}
	if (HasStop)
		{
		uint StopAa = AaStart + AaLength;
		uint Pos = T.GetNtPos(StopAa, 1);
		asserta(Pos < NtL);
		asserta(Row[Pos] == ' ');
		Row[Pos] = (char) T.GetAa(StopAa);
		}
	return Row;
	}

static void AppendFrameOrfs(const Translator &T, const string &Label,
  const string &Sep, uint &OrfIndex, bool EndlessOnly,
  bool PlusOnly, bool MinusOnly, bool WantAa, bool WantNt, bool WantAln,
  vector<OrfHit> &Hits)
	{
	const uint NtL = T.m_NtL;
	const vector<uint> &AaStarts = T.m_OrfAaStarts;
	const vector<uint> &AaLengths = T.m_OrfAaLengths;
	const vector<bool> &HasStops = T.m_OrfHasStops;
	const vector<bool> &HasMets = T.m_OrfHasMets;
	const vector<bool> &EndlessVec = T.m_OrfEndless;

	const uint N = SIZE(AaStarts);
	asserta(SIZE(AaLengths) == N);
	asserta(SIZE(HasStops) == N);
	asserta(SIZE(HasMets) == N);
	asserta(SIZE(EndlessVec) == N);

	int Frame = T.m_Frame;
	if (PlusOnly && Frame < 0)
		return;
	if (MinusOnly && Frame > 0)
		return;

	for (uint i = 0; i < N; ++i)
		{
		uint AaStart = AaStarts[i];
		uint AaLength = AaLengths[i];
		bool HasStop = HasStops[i];
		bool HasMet = HasMets[i];
		bool Endless = EndlessVec[i];

		if (EndlessOnly && !Endless)
			continue;

		uint NtStart = T.GetNtPos(AaStart, 0);
		uint NtEnd = NtStart + 3*AaLength - 1;

		++OrfIndex;
		OrfHit Hit;
		Hit.AaLength = AaLength;
		Hit.Frame = Frame;
		Hit.Label = Label;
		Hit.Label += Sep;
		Psa(Hit.Label, "orf%u", OrfIndex);
		Hit.Label += Sep;
		Psa(Hit.Label, "code%u", GetNCBICodeNr());
		Hit.Label += Sep;
		Psa(Hit.Label, "%u-%u/%u(%+d)",
		  NtStart + 1, NtEnd + 1, NtL, Frame);
		if (Endless)
			{
			if (NtL%3 == 0)
				Hit.Label += "Ce";
			else
				Hit.Label += "Co";
			if (!HasMet)
				Hit.Label += "m";
			}

		if (WantAa)
			T.GetAaSubStr(AaStart, AaLength, Hit.AaSeq);
		if (WantNt)
			T.GetNtSubStr(AaStart, AaLength, Hit.NtSeq);
		if (WantAln)
			Hit.AlnRow = MakeAlnRow(T, AaStart, AaLength, HasStop);
		Hits.push_back(Hit);
		}
	}

// Longer codon count wins. Equal length prefers plus strand over minus.
// Same strand keeps the earlier hit.
static void KeepLongestOrf(vector<OrfHit> &Hits)
	{
	const uint N = SIZE(Hits);
	if (N <= 1)
		return;

	uint Best = 0;
	for (uint i = 1; i < N; ++i)
		{
		const OrfHit &Cand = Hits[i];
		const OrfHit &Win = Hits[Best];
		if (Cand.AaLength > Win.AaLength)
			Best = i;
		else if (Cand.AaLength == Win.AaLength &&
		  Cand.Frame > 0 && Win.Frame < 0)
			Best = i;
		}
	if (Best != 0)
		Hits[0] = Hits[Best];
	Hits.resize(1);
	}

static void WriteAlnout(FILE *f, const byte *Nt, uint NtL, const string &Label,
  const vector<OrfHit> &Hits, uint Cols)
	{
	if (f == 0 || NtL == 0)
		return;

	uint Block = Cols;
	if (Block == 0)
		Block = NtL;

	for (uint i = 0; i < SIZE(Hits); ++i)
		asserta(SIZE(Hits[i].AlnRow) == NtL);

	for (uint From = 0; From < NtL; )
		{
		uint n = NtL - From;
		if (n > Block)
			n = Block;
		fprintf(f, "%.*s  %s\n", (int) n, (const char *) (Nt + From),
		  Label.c_str());
		for (uint i = 0; i < SIZE(Hits); ++i)
			fprintf(f, "%.*s  %s\n", (int) n,
			  Hits[i].AlnRow.c_str() + From, Hits[i].Label.c_str());
		fprintf(f, "\n");
		From += n;
		}
	}

void RunTranslate(const string &InputFileName)
	{
	const string &OutputFileName = opt(output);

	if (optset_aaout)
		Die("Use -output for aa seqs");

	FILE *fOut = CreateStdioFile(OutputFileName);
	FILE *fNtOut = CreateStdioFile(opt(ntout));

	FASTASeqSource SS;
	SS.Open(InputFileName);

	bool Circular = opt(circular);
	if (optset_mincodons)
		opt(mincodons);

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
			string AaStr;
			T.GetAaStr(AaStr);

			string OutLabel;
			OutLabel = Label;
			Psa(OutLabel, " frame%+d", Frame);
			SeqToFasta(fOut, (const byte *) AaStr.c_str(), SIZE(AaStr),
			  OutLabel.c_str());
			}
		}

	CloseStdioFile(fOut);
	CloseStdioFile(fNtOut);
	}

void RunGetOrfs(const string &InputFileName, bool DefaultCircular)
	{
	if (optset_aaout)
		Die("Use -output for aa seqs");

	bool SingleEndless = opt(single_endless);
	bool DoubleEndless = opt(double_endless);
	asserta(!(SingleEndless && DoubleEndless));

	if (opt(circular) && opt(linear))
		Die("-circular and -linear are incompatible");
	if ((SingleEndless || DoubleEndless) && opt(linear))
		Die("-linear cannot be combined with -single_endless or -double_endless");

	bool Circular = DefaultCircular;
	if (optset_circular)
		Circular = true;
	if (opt(linear))
		Circular = false;
	if (SingleEndless || DoubleEndless)
		Circular = true;

	bool PlusOnly = opt(plus_only);
	bool MinusOnly = opt(minus_only);
	bool LongestOrf = opt(longest_orf);
	uint MinCodonCount = opt(mincodons);

	string Sep = " ";
	if (optset_sep)
		Sep = opt(sep);

	uint AlnCols = 100;
	if (optset_fasta_cols)
		AlnCols = opt(fasta_cols);

	FILE *fAa = CreateStdioFile(opt(output));
	FILE *fNt = CreateStdioFile(opt(ntout));
	FILE *fFev = CreateStdioFile(opt(fevout));
	FILE *fAln = CreateStdioFile(opt(alnout));

	FASTASeqSource SS;
	SS.Open(InputFileName);
	uint64 FileSize = GetStdioFileSize64(SS.m_LR.m_f);

	vector<int> Frames;
	GetFrames(Frames);

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	Translator T;

	bool WantAa = (fAa != 0);
	bool WantNt = (fNt != 0);
	bool WantAln = (fAln != 0);
	bool EndlessOnly = (SingleEndless || DoubleEndless);

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

		const byte *Q = SIQ->m_Seq;
		const uint QL = SIQ->m_L;
		const string Label = string(SIQ->m_Label);

		T.SetSeq(Q, QL, Circular);

		if (EndlessOnly)
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

		if (fFev != 0 && !LongestOrf)
			{
			string FevStr;
			GetFevStr(SIQ, Circular, Frames, MinCodonCount,
			  PlusOnly, MinusOnly, FevStr);
			if (FevStr != "")
				fprintf(fFev, "%s\t%s\n", SIQ->m_Label, FevStr.c_str());
			}

		if (WantAa || WantNt || WantAln || LongestOrf)
			{
			vector<OrfHit> Hits;
			uint OrfIndex = 0;
			for (uint k = 0; k < SIZE(Frames); ++k)
				{
				int Frame = Frames[k];
				if (!T.ValidFrame(Frame))
					continue;

				T.SetFrame(Frame);
				T.SetOrfs(MinCodonCount);
				AppendFrameOrfs(T, Label, Sep, OrfIndex, EndlessOnly,
				  PlusOnly, MinusOnly, WantAa, WantNt, WantAln, Hits);
				}
			if (LongestOrf)
				KeepLongestOrf(Hits);
			for (uint i = 0; i < SIZE(Hits); ++i)
				{
				if (WantAa)
					SeqToFasta(fAa, Hits[i].AaSeq.c_str(), SIZE(Hits[i].AaSeq),
					  Hits[i].Label.c_str());
				if (WantNt)
					SeqToFasta(fNt, Hits[i].NtSeq.c_str(), SIZE(Hits[i].NtSeq),
					  Hits[i].Label.c_str());
				}
			if (WantAln)
				WriteAlnout(fAln, Q, QL, Label, Hits, AlnCols);
			if (LongestOrf && fFev != 0 && SIZE(Hits) > 0)
				{
				uint OrfLength = Hits[0].AaLength;
				double Pct = (100.0*OrfLength)/double(QL);
				string FevStr;
				Ps(FevStr, "maxorf=%u/%.1f", OrfLength, Pct);
				fprintf(fFev, "%s\t%s\n", SIQ->m_Label, FevStr.c_str());
				}
			}
		}
	ProgressStep(999, 1000, "Processing (%s seqs)", IntToStr(SeqCount));

	CloseStdioFile(fAa);
	CloseStdioFile(fNt);
	CloseStdioFile(fFev);
	CloseStdioFile(fAln);
	}

void cmd_translate()
	{
	RunTranslate(opt(translate));
	}

void cmd_getorfs()
	{
	RunGetOrfs(opt(getorfs), false);
	}
