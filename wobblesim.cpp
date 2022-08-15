#include "myutils.h"
#include "translator.h"
#include "seqinfo.h"
#include "objmgr.h"

void GetRandEndlessOrf(byte *Seq, uint CodonCount);
void CountWobblemers(SeqInfo *SI, vector<uint> &Counts);
void WriteWobblemerCounts(FILE *f, const vector<uint> &Counts);

void cmd_wobblesim()
	{
	FILE *fOut = CreateStdioFile(opt(wobblesim));
	const uint IterCount = 10000000;

	const uint DictSize = 1024; // myipow(4, 5);
	vector<uint> Counts(DictSize, 0);

	Translator T;

	const uint MinCodonCount = 75;
	const uint MaxCodonCount = 150;
	const uint d = MaxCodonCount - MinCodonCount;
	byte *Seq = myalloc(byte, MaxCodonCount*3);
	vector<uint> StopCounts(3, 0);
	vector<uint> LetterCounts(3, 0);
	uint DoubleCount = 0;
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	for (uint Iter = 0; Iter < IterCount; ++Iter)
		{
		ProgressStep(Iter, IterCount, "Simulating, %u doubles", DoubleCount);
		uint CodonCount = MinCodonCount + randu32()%d;
		asserta(CodonCount >= MinCodonCount);
		asserta(CodonCount <= MaxCodonCount);
		GetRandEndlessOrf(Seq, CodonCount);

		uint L = CodonCount*3;
		T.SetSeq(Seq, L, true);

		int PF = 0;
		int MF = 0;

		T.GetEndlessFrames2(PF, MF);
		asserta(PF != 0);
		if (MF == 0)
			continue;
		++DoubleCount;

		SI->SetPtrs(0, "rand", Seq, L);
		CountWobblemers(SI, Counts);
		}

	ProgressLog("%u doubles\n", DoubleCount);
	WriteWobblemerCounts(fOut, Counts);
	}
