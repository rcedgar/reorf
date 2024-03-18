#include "myutils.h"
#include "seqdb.h"
#include "translator.h"
#include <map>

static byte GetMut(byte c)
	{
	for (uint i = 0; i < 1000; ++i)
		{
		uint k = randu32()%4;
		byte Mutc = "ACGT"[k];
		if (Mutc != c)
			return Mutc;
		}
	asserta(false);
	return 0;
	}

static bool DoIter(const byte *Seq, uint L, uint PctId)
	{
	byte *MutSeq = myalloc(byte, L);
	memcpy(MutSeq, Seq, L);

	uint MutCount = ((100 - PctId)*L)/100;
	for (uint MutIndex = 0; MutIndex < MutCount; ++MutIndex)
		{
		uint Pos = randu32()%L;
		uint k = randu32()%3;
		char c = MutSeq[Pos];
		byte Mutc = GetMut(c);
		MutSeq[Pos] = Mutc;
		}

	Translator T;
	T.SetSeq(MutSeq, L, true);
	bool MutHasStop = false;
	if (opt(double_endless))
		{
		int PlusFrame, MinusFrame;
		T.GetEndlessFrames2(PlusFrame, MinusFrame);
		MutHasStop = !(PlusFrame > 0 && MinusFrame < 0);
		}
	else
		MutHasStop = !T.CurrentFrameIsEndless();

	string MutAa;
	T.GetAaStr(MutAa);
	myfree(MutSeq);
	return MutHasStop;
	}

static double DoSeq(const byte *Seq, uint L, uint PctId, uint Iters)
	{
	uint StopCount = 0;
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		bool HasStop = DoIter(Seq, L, PctId);
		if (HasStop)
			++StopCount;
		}
	return double(StopCount)/Iters;
	}

void cmd_evolve_neutral()
	{
	SeqDB Input;
	Input.FromFasta(opt(evolve_neutral));
	const uint SeqCount = Input.GetSeqCount();

	uint Iters = 10000;
	if (optset_iters)
		Iters = opt(iters);
	ProgressLog("Iters %u\n", Iters);

	vector<uint> PctIds;
	PctIds.push_back(99);
	PctIds.push_back(98);
	PctIds.push_back(95);
	PctIds.push_back(90);
	PctIds.push_back(80);
	PctIds.push_back(70);
	const uint PctIdCount = SIZE(PctIds);

	vector<uint> Ls;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const uint L = Input.GetSeqLength(SeqIndex);
		Ls.push_back(L);
		}

	Translator T;

	map<pair<uint, uint>, double> KeyToProb;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Evolving");

		const byte *Seq = Input.GetSeq(SeqIndex);
		const uint L = Input.GetSeqLength(SeqIndex);
		T.SetSeq(Seq, L, true);
		T.SetFrame(1);
		asserta(T.CurrentFrameIsEndless());
		for (uint PctIdIndex = 0; PctIdIndex < PctIdCount; ++PctIdIndex)
			{
			uint PctId = PctIds[PctIdIndex];
			double P_stop = DoSeq(Seq, L, PctId, Iters);
			pair<uint, uint> Key(L, PctId);
			KeyToProb[Key] = P_stop;
//			Log("Seq %u PctId = %u L = %u P_stop = %.4f\n", SeqIndex, PctId, L, P_stop);
			}
		}

	Log("\n\n");
	Log("P_stop");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint L = Ls[SeqIndex];
		string s;
		Ps(s, "L=%u", L);
		Log("  %7.7s", s.c_str());
		}
	Log("\n");

	for (uint PctIdIndex = 0; PctIdIndex < PctIdCount; ++PctIdIndex)
		{
		uint PctId = PctIds[PctIdIndex];
		Log("%5u%%", PctId);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			uint L = Ls[SeqIndex];
			pair<uint, uint> Key(L, PctId);
			double P_stop = KeyToProb[Key];
			Log("   %6.3f", P_stop);
			}
		Log("\n");
		}
	}
