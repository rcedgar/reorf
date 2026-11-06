#include "myutils.h"
#include "seqinfo.h"
#include "translator.h"

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

static bool FrameSelected(int Frame, const vector<int> &Frames,
  bool PlusOnly, bool MinusOnly)
	{
	if (PlusOnly && Frame < 0)
		return false;
	if (MinusOnly && Frame > 0)
		return false;
	for (uint i = 0; i < SIZE(Frames); ++i)
		{
		if (Frames[i] == Frame)
			return true;
		}
	return false;
	}

void GetFevStr(const SeqInfo *SIQ, bool Circular, const vector<int> &Frames,
  uint MinCodonCount, bool PlusOnly, bool MinusOnly, string &Str)
	{
	Str.clear();

	const byte *Q = SIQ->m_Seq;
	const uint QL = SIQ->m_L;

	Translator T;
	T.SetSeq(Q, QL, Circular);

	vector<int> EndlessFrames;
	vector<bool> EndlessHasMets;
	T.GetEndlessFrames(EndlessFrames, EndlessHasMets);

	vector<int> KeptFrames;
	vector<bool> KeptMets;
	for (uint i = 0; i < SIZE(EndlessFrames); ++i)
		{
		int Frame = EndlessFrames[i];
		if (!FrameSelected(Frame, Frames, PlusOnly, MinusOnly))
			continue;
		KeptFrames.push_back(Frame);
		KeptMets.push_back(EndlessHasMets[i]);
		}

	const uint F = SIZE(KeptFrames);
	if (F == 0)
		{
		bool Unrestricted = (!PlusOnly && !MinusOnly && !optset_frame);
		if (Unrestricted)
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

		uint MaxLen = 0;
		bool Any = false;
		for (uint k = 0; k < SIZE(Frames); ++k)
			{
			int Frame = Frames[k];
			if (!T.ValidFrame(Frame))
				continue;
			if (!FrameSelected(Frame, Frames, PlusOnly, MinusOnly))
				continue;
			T.SetFrame(Frame);
			T.SetOrfs(MinCodonCount);
			for (uint j = 0; j < SIZE(T.m_OrfAaLengths); ++j)
				{
				uint Len = T.m_OrfAaLengths[j];
				if (Len > MaxLen)
					{
					MaxLen = Len;
					Any = true;
					}
				}
			}
		if (Any)
			{
			double Pct = (100.0*MaxLen)/double(QL);
			Ps(Str, "maxorf=%u/%.1f", MaxLen, Pct);
			}
		return;
		}

	uint NEpYesMet = 0;
	uint NEmYesMet = 0;
	uint NEpNoMet = 0;
	uint NEmNoMet = 0;
	for (uint i = 0; i < F; ++i)
		{
		int Frame = KeptFrames[i];
		bool HasMet = KeptMets[i];
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

void RunGetOrfs(const string &InputFileName, bool DefaultCircular);

void cmd_getorfs_serratus()
	{
	RunGetOrfs(opt(getorfs_serratus), true);
	}
